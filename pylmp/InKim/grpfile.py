#!/home/noische/Enthought/Canopy_64bit/User/bin/python
import nutils as nu

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

class grpfile():

    #grp[group_no]: atoms, natoms, constraints, rotsym, linear, energy, volume

    def __init__(self, grp_file):
        print("Reading group file %s" % grp_file)

        #self.grp = nu.hash()
        self.grp = AutoVivification()
        with open(grp_file) as f:
            for line in f:
                # Total group numbers
                if "Total Groups:" in line:
                    n_group_in_file = int(line.strip("Total Groups:"))
                elif "Group " in line:
                    # group definition
                    parse = line.split()
                    try:
                        group_no = int(parse[1])
                        group_natoms = int(parse[3])
                    except ValueError:
                        nu.die("Error on reading group definitions!")

                    self.grp[group_no]['natoms'] = group_natoms # number of atoms declared

                    # group member atoms
                    atomline = next(f)
                    atomline = atomline.replace(' - ', '-')
                    atomline = atomline.replace(' ', '')
                    self.grp[group_no]['atoms'] = nu.range_expand(atomline) # atoms specified

                elif "Constraints" in line:
                    parse = next(f).split()
                    for index, i in enumerate(parse):
                        try:
                            self.grp[index+1]['constraints'] = int(i)
                        except ValueError:
                            nu.die("Error on reading constraints %s!" % i)

                elif "RotationalSymmetryNumber" in line:
                    parse = next(f).split()
                    for index, i in enumerate(parse):
                        try:
                            self.grp[index+1]['rotsym'] = int(i)
                        except ValueError:
                            nu.die("Error on reading rotational symmetry %s!" % i)
                elif "LinearMoleculeFlag" in line:
                    parse = next(f).split()
                    for index, i in enumerate(parse):
                        try:
                            self.grp[index+1]['linear'] = int(i)
                        except ValueError:
                            nu.die("Error on reading linear molecular flag %s!" % i)
                elif "GroupEnergyAvg" in line:
                    parse = next(f).split()
                    for index, i in enumerate(parse):
                        try:
                            self.grp[index+1]['energy'] = float(i)
                        except ValueError:
                            nu.die("Error on reading linear molecular flag %s!" % i)
                elif "GroupVolume" in line:
                    parse = next(f).split()
                    for index, i in enumerate(parse):
                        try:
                            self.grp[index+1]['volume'] = float(i)
                        except ValueError:
                            nu.die("Error on reading linear molecular flag %s!" % i)

        # integrity check: group definition mismatch
        if n_group_in_file != len(self.grp.keys()):
            nu.die("Number of groups mismatch! %d total groups declared but %d groups defined." % (n_group_in_file, len(self.grp.keys())))
        # integrity check: wrong group numbering
        if len(self.grp.keys()) != self.grp.keys()[-1]:
            nu.die("Number of groups mismatch! %d total groups declared but %d groups not specified." % (n_group_in_file, len(self.grp.keys()) - self.grp.keys()[-1]))
        # integrity check: wrong atom number specified
        for gno in self.grp.keys():
            if len(self.grp[gno]['atoms']) != self.grp[gno]['natoms']:
                nu.die("Declared atoms and specified atom numbers mismatch in group %d." % group_no)

        print("Finished reading group file %s" % grp_file)


    def write(self, out_file, zip=False):
        if zip:
            nu.warn("Atoms will be specified in range if possible. (i.e. xx - yy)")
        else:
            nu.warn("Atom numbers will be recorded in full specification.")

        with open(out_file, 'w') as f:
            n_group = len(self.grp.keys())
            groups = self.grp.keys()
            groups.sort()
            
            f.write("Total Groups: %d\n" % n_group)
            
            for i in groups:
                atoms = sorted(self.grp[i]['atoms'])
                f.write("Group %d Atoms %d\n" % (i, len(atoms)))
                if zip:
                    f.write(''.join((('%i - %i ' % r) if len(r) == 2 else '%i' % r) for r in nu.range_extract(atoms)) + "\n")
                else:
                    f.write(''.join("%d " % i for i in atoms) + "\n")

            f.write("Constraints" + "\n")
            f.write(''.join("%d " % self.grp[i]['constraints'] for i in groups) + "\n")

            f.write("RotationalSymmetryNumber" + "\n")
            f.write(''.join("%d " % self.grp[i]['rotsym'] for i in groups) + "\n")

            f.write("LinearMoleculeFlag" + "\n")
            f.write(''.join("%d " % self.grp[i]['linear'] for i in groups) + "\n")

            if self.grp[1]['volume']:
                f.write("GroupVolume" + "\n")
                f.write(''.join("%-16.5f " % self.grp[i]['volume'] for i in groups) + "\n")

            if self.grp[1]['energy']:
                f.write("GroupEnergy" + "\n")
                f.write(''.join("%-16.5f " % self.grp[i]['energy'] for i in groups) + "\n")

    def split_group(self, group_no, lst):
        """
        Splits a new group with specified atoms in the list lst from self.grp[group_no]
        """
        # returns nothing with an empty list
        if not lst:
            print("The group cannot be separated into an empty group!")
            return -1

        # check if atoms in lst also exists in group_no
        for i in lst:
            if not i in self.grp[group_no]['atoms']:
                print("Atom %d does not exist in the group %d! The group cannot be split." % (i, group_no))
                return -1

        # if specified lst is the whole atoms list, no need to split the group
        if set(self.grp[group_no]['atoms']) == set(lst) or len(lst) == len(self.grp[group_no]['atoms']):
            print("The group cannot be separated with same atom numbers!")
            return -1

        # consider constraints
        original_n_constraints = self.grp[group_no]['constraints']
        original_n_atoms = len(self.grp[group_no]['atoms'])
        q = original_n_constraints / original_n_atoms

        subs = [i for i in self.grp[group_no]['atoms'] if not i in lst]
        max_group_no = max(self.grp.keys())
        new_group_no = max_group_no + 1
        
        # copy
        self.grp[group_no]['atoms'] = subs
        self.grp[new_group_no]['atoms'] = lst

        if q: 
            self.grp[group_no]['constraints'] = len(self.grp[group_no]['atoms']) / q
            self.grp[new_group_no]['constraints'] = len(self.grp[new_group_no]['atoms']) / q

        self.grp[new_group_no]['rotsym'] = self.grp[group_no]['rotsym']
        self.grp[new_group_no]['linear'] = self.grp[group_no]['linear']
        if self.grp[group_no]['volume']: self.grp[new_group_no]['volume'] = self.grp[group_no]['volume']
        if self.grp[group_no]['energy']: self.grp[new_group_no]['energy'] = self.grp[group_no]['energy']

        return new_group_no

    def find_group(self, ano):
        """
        Returns group numbers which contains ano.
        """
        result = [];
        for key in self.grp.keys():
            if ano in self.grp[key]['atoms']:
                result.append(key)
        
        return result

    def copy_group(self, group_no):
        """
        Duplicate a new group.
        """
        max_group_no = max(self.grp.keys())
        new_group_no = max_group_no + 1

        self.grp[new_group_no]['atoms'] = self.grp[group_no]['atoms']
        self.grp[new_group_no]['constraints'] = self.grp[group_no]['constraints']
        self.grp[new_group_no]['rotsym'] = self.grp[group_no]['rotsym']
        self.grp[new_group_no]['linear'] = self.grp[group_no]['linear']
        if self.grp[group_no]['volume']: self.grp[new_group_no]['volume'] = self.grp[group_no]['volume']
        if self.grp[group_no]['energy']: self.grp[new_group_no]['energy'] = self.grp[group_no]['energy']

        return new_group_no
