#!/home/noische/python
"""
CNT.py

Module containing CNT-related tools including:
(*) CNT class

20160120 Now we use atom.chain to determine internal and exterior water molecules.
"""

import os, sys, copy
import numpy as np
import bgf, bgftools
import nutils as nu
import random

version = "150204"

class Nanotube(object):

    bgfmodel = bgf.BgfFile()
    type = None;
    n_allatoms = 0;
    n_nanotube = 0;
    n_water = 0;
    NTatoms = [];
    WATatoms = [];
    otheratoms = [];

    isAligned = None;
    isPBCadjusted = None;

    isRadiusCalc = None;
    orientation = None;    # nanotube direction
    axis = None;    # nparray mask according to nanotube direction. x direction = [0, 1, 1] for numpy
    axismask = None;    # nparray mask for nanotube direction. x direction = [1, 0, 0]
    radius = None;
    height = None;
    zhi = None;
    zlo = None;
    water_position_determined = None;    # True if determine_water_position() performed
    n_inside_water = None;
    n_outside_water = None;

    pbc = [];   # [xhi, yhi, zhi]
    ff = "";
    model_filename = "";

    def __init__(self, fname, *args, **kwargs):

        bgf_file = "";
        if fname:
            bgf_file = fname

            if isinstance(bgf_file, bgf.BgfFile):
                self.bgfmodel = bgf_file
                self.model_filename = ""    # if the model is directly from bgf object, no filename will be assigned.
            else:
                self.bgfmodel = bgf.BgfFile(bgf_file);
                self.model_filename = bgf_file

            #
            # Determine nanotube type and atom types
            #
            self.NTatoms = self.find_NT_atoms();
            self.WATatoms = self.find_WAT_atoms();
            self.otheratoms = self.find_other_atoms();
            self.find_type();
            self.set_ff();
            
            #
            # Set number of atoms
            # 
            self.n_allatoms = len(self.bgfmodel.a)
            self.n_nanotube = len(self.NTatoms)
            self.n_water = len(self.WATatoms)

            #
            # Properties of CNT: requires special treatment so not automatically calculated.
            #
            self.isAligned = False;
            self.isRadiusCalc = False;
            self.isPBCadjusted = False;
            self.orientation = None;
            self.radius = None;
            self.height = None;
            self.zhi = None;
            self.zlo = None;
            self.water_position_determined = False;

            # model check-up: total charge
            chg = bgftools.charge(self.bgfmodel)
            if chg != 0 and abs(chg) > 1e-10:
                nu.warn("Charge is not neutral: " + str(chg))

            if self.bgfmodel.CRYSTX != []:
                self.pbc = self.bgfmodel.CRYSTX[:3]
        else:
            nu.warn("Nanotube class: An empty object created.")

    def update(self):
        self.NTatoms = self.find_NT_atoms();
        self.WATatoms = self.find_WAT_atoms();
        self.otheratoms = self.find_other_atoms();
        self.n_allatoms = len(self.bgfmodel.a)
        self.n_nanotube = len(self.NTatoms)
        self.n_water = len(self.WATatoms)


    #
    # nanotube atoms
    #
    def find_NT_atoms(self):
        _NT_atoms = [];
        for atom in self.bgfmodel.a:
            if "NT" in atom.rName:
                _NT_atoms.append(atom)
        if _NT_atoms == []:
            nu.warn("find_NT_atoms: No atoms with residue name NT found.")

        print("find_NT_atoms(): found " + str(len(_NT_atoms)) + " atoms for the nanotube in BGF file.")
        return _NT_atoms;


    def find_WAT_atoms(self):
        _WAT_atoms = [];
        for atom in self.bgfmodel.a:
            if "WAT" in atom.rName:
                _WAT_atoms.append(atom)
        if _WAT_atoms == []:
            nu.warn("find_WAT_atoms: No atoms with residue name WAT found.")

        print("find_WAT_atoms(): found " + str(len(_WAT_atoms)) + " atoms for water in BGF file. (" + str(len(_WAT_atoms)/3) + " molecules)")
        return _WAT_atoms;

    def find_other_atoms(self):
        _other_atoms = []
        for atom in self.bgfmodel.a:
            if "NT" in atom.rName or "WAT" in atom.rName:
                pass;
            else:
                _other_atoms.append(atom)

        print("find_other_atoms(): found " + str(len(_other_atoms)) + " atoms for other atoms in BGF file." )
        return _other_atoms;

    #
    # nanotube type
    #
    def set_type(self, _type=""):
        if _type != "":
            nu.warn("set_type(): You are manually assigning the nanotube type.")
            self.type = _type


    def get_type(self):
        return self._type


    def find_type(self):
        #_NT_atoms = self.find_NT_atoms()
        if len(self.NTatoms) != 0:
            self.type = self.NTatoms[0].rName
        else:
            nu.warn("find_type(): Failed to find Nanotube types.")

    ### end of nanotube type ###


    #
    # nanotube BGF model
    #
    def set_ff(self, ff_file=""):
        if ff_file != "":
            self.ff = ff_file
        if not self.type:
            self.find_type()

        if self.type == "CNT":
            self.ff = "/home/noische/ff/graphite_lj.ff"
        if self.type == "BNT":
            self.ff = "/home/noische/ff/DREIDING2.21.ff"


    def update_pbc(self, pbc):
        if len(pbc) < 3:
            nu.warn("update_pbc(): Wrong assignment for the dimension of the model. Must be 3.")
        elif len(pbc) == 3:
            self.bgfmodel.CRYSTX = pbc + [90.0, 90.0, 90.0]
            self.pbc = pbc
        elif len(pbc) == 6:
            self.bgfmodel.CRYSTX = pbc
            self.pbc = pbc[:3]


    ### end of nanotube orientation ###

                
    #
    # functions
    #
    def coord2nparray(self, bgfatoms):
        temp = [];
        for atom in bgfatoms:
            temp.append([atom.x, atom.y, atom.z])
        return np.array(temp)


    def set_coord(self, coord):
        if len(coord) != len(self.bgfmodel.a):
            nu.die("Number of coordinates != number of BGF atoms")

        for index, i in enumerate(coord):
            self.bgfmodel.a[index].x = i[0]
            self.bgfmodel.a[index].y = i[1]
            self.bgfmodel.a[index].z = i[2]

        
    def make_infinite(self, *args):
        '''
        args: filename to save
        '''
        print("CNT.py: make_infinite()")
        if len(args) > 1:
            nu.warn("make_infinite(): Too many arguments.")
            
        if not self.isPBCadjusted: 
            nu.warn("make_infinite(): You are trying to make infinite NT without adjusting PBC.")
            #return 0;

        print("\t* Found " + str(len(self.NTatoms)) + " atoms for the nanotube.")
        print("\t* Found " + str(len(self.WATatoms)) + " atoms for water.")


        if self.type == "BNT":
            self.detach_hydrogen()

        aNo_pair = [];
        for atom in self.NTatoms:
            if len(atom.CONECT) == 3:
                continue;

            x = np.array([atom.x, atom.y, atom.z])
            min_atom_aNo = 100000;    # placeholder
            min_atom_d = 10000.0;

            for atom2 in self.NTatoms:
                if (not atom2.aNo in atom.CONECT) and len(atom2.CONECT) == 2:
                    if (atom.rName == atom2.rName) and (self.type == "CNT" or atom.ffType != atom2.ffType):
                        y = np.array([atom2.x, atom2.y, atom2.z])
                        if -1.0 < ((x - y) * self.axismask).sum() < 1.0:
                            continue;    # prevent adjacent atoms
                        d = nu.pbc_dist(x, y, self.pbc)
                        if d < min_atom_d:
                            min_atom_aNo = atom2.aNo
                            min_atom_d = d
                        
            aNo_pair.append([atom.aNo, min_atom_aNo])

        for i in aNo_pair:
            self.bgfmodel.connectAtoms(i[0], i[1])

        self.check_connectivity();
        self.bgfmodel.renumber()

        if len(args) == 0:
            value = raw_input("Do you want to save the infinite Nanotube structure to BGF file [N]? ")
            if "y" in value.lower():
                filename = self.model_filename[:-4] + ".infinite.bgf"
                value = raw_input("Filename to save [" + filename + "]? ") or filename
                self.bgfmodel.saveBGF(value)
        elif len(args) == 1:
            self.bgfmodel.saveBGF(args[0])


        ### end of make_infinite


    def check_connectivity(self):
        # connectivity check-up
        print("Checking nanotube connectivity..")
        try:
            for atom in self.NTatoms:
                if len(atom.CONECT) != 3:
                    raise ValueError
        except ValueError:
            nu.warn("Defect on nanotube connection found: " + str(atom.CONECTline()))
        else:
            print("  No connectivity error found.")
            

    def calc_height_radius(self):
        """
        sets center, height, orientation, axis
        """
        if not self.isAligned:
            nu.warn("calc_height_radius: Nanotube is not aligned along axes. Values may be inaccurate.")

        nt_coord = self.coord2nparray(self.NTatoms)

        min_x, min_y, min_z = nt_coord.min(axis=0)
        max_x, max_y, max_z = nt_coord.max(axis=0)

        range = np.array([[min_x, max_x], [min_y, max_y], [min_z, max_z]])
        diff = [max_x - min_x, max_y - min_y, max_z - min_z]

        #self.center = np.mean(nt_coord, axis=0)        # set nanotube center
        self.center = self.get_nt_center();        # set nanotube center
        #std_c = np.std(nt_coord, axis=0)

        # height and orientation estimation: can be inaccurate
        self.height = np.max(diff)        # set nanotube height
        if np.argmax(diff) == 0:
            self.orientation = "x"
            self.axis = np.array([0, 1, 1])
            self.axismask = np.array([1, 0, 0])
        if np.argmax(diff) == 1:
            self.orientation = "y"
            self.axis = np.array([1, 0, 1])
            self.axismask = np.array([0, 1, 0])
        if np.argmax(diff) == 2:
            self.orientation = "z"
            self.axis = np.array([1, 1, 0])
            self.axismask = np.array([0, 0, 1])

        self.zlo, self.zhi = np.sum(range.T * self.axismask, axis=1)
        self.radius = np.mean(np.sqrt(((nt_coord * self.axis - self.center * self.axis)**2).sum(axis=-1)))    # set nanotube radius
        self.bgfmodel.REMARK.append("Nanotube center: %s" % self.center)
        self.bgfmodel.REMARK.append("Nanotube height: %8.3f" % self.height)
        self.bgfmodel.REMARK.append("Nanotube radius: %8.3f" % self.radius)

        print("\t* Nanotube height: "+ str(self.height))
        print("\t* Nanotube radius: "+ str(self.radius))

        self.isRadiusCalc = True


    def get_nt_center(self):
        nt_coord = self.coord2nparray(self.NTatoms)
        return np.mean(nt_coord, axis=0)


    def adjust_pbc(self, out_file=""):
        """
        adjust pbc of BGF model for infinite connectivity.
        """
        if not self.isRadiusCalc:
            nu.warn("adjust_pbc: You are adjusting PBC without calculating height and radius of CNT.")
            self.calc_height_radius()

        # set PBC
        _dim = self.pbc
        
        margin = 0.0
        val = raw_input("Do you want to change the margin [N]? ")
        if "y" in val.lower():
            val = raw_input("Margin from the radius to the pbc [%3.1f]? " % margin) or 0.0

            if self.orientation == "x":
                _dim = [self.height + 1.0, self.radius + margin, self.radius + margin]
            elif self.orientation == "y":
                _dim = [self.radius + margin, self.height + 1.0, self.radius + margin]
            elif self.orientation == "z":
                _dim = [self.radius + margin, self.radius + margin, self.height + 1.0]

        else:
            if self.orientation == "x":
                _dim[0] = self.height + 1.0
            if self.orientation == "y":
                _dim[1] = self.height + 1.0
            if self.orientation == "z":
                _dim[2] = self.height + 1.0

        self.update_pbc(_dim)

        # move CM to boxcenter
        boxcenter = np.array(self.pbc) / 2
        cm = np.array(self.get_nt_center())
        delta = boxcenter - cm
        coord = self.coord2nparray(self.bgfmodel.a)

        self.set_coord(coord + delta)

        self.isPBCadjusted = True


    def make_pbc(self, *args):
        """
        assign pbc information on BGF model.
        """
        self.bgfmodel.PERIOD = "111"
        self.bgfmodel.AXES = "ZYX"
        self.bgfmodel.SGNAME = "P 1                  1    1"
        self.bgfmodel.CELLS = [-1, 1, -1, 1, -1, 1]

        if len(args) == 1:
            if len(args[0]) == 6:
                self.bgfmodel.CRYSTX = args[0]
            elif len(args[0]) == 3:
                self.bgfmodel.CRYSTX = args[0] + [90.0, 90.0, 90.0]

        if self.bgfmodel.CRYSTX == []:
            self.adjust_pbc()


    def make_bulk_infinite(self, *args):
        """
        calculate nanotube height, trim water molecules outside the box,
        and connect the two ends to make infinite nanotube.
        """
        self.make_pbc()
        self.calc_height_radius()
        self.adjust_pbc()

        # trim
        temp_list = []
        for atom in self.bgfmodel.a:
            if "O" in atom.ffType and "WAT" in atom.rName:
                if atom.z > self.pbc[2] or atom.z < 0.0:
                    temp_list.append(atom.aNo)
                    temp_list += atom.CONECT

        temp_list = list(set(temp_list))
        print("\tFound " + str(len(temp_list)) + " atoms for water outside of the box: will be removed.")
        del_list = []
        for i in temp_list:
            del_list.append(self.bgfmodel.a2i[i])
        self.bgfmodel.delAtoms(del_list, False);
        self.bgfmodel.renumber();
        print('\tTrim successful!')
        self.update()

        # connect
        if len(args) == 0:
            self.make_infinite()
        elif len(args) == 1:
            self.make_infinite(args[0])
        else:
            nu.warn("Wrong arguments passed.")
            return 0;


    def remove_exterior_water(self, *args):
        temp_list = [];
        if not self.water_position_determined:
            self.determine_water_position()

        # delete WTO atoms
        for atom in self.bgfmodel.a:
            if "O" in atom.chain and "WAT" in atom.rName:
                temp_list.append(atom.aNo)
                temp_list += atom.CONECT
        temp_list = list(set(temp_list))
        del_list = []
        for i in temp_list:
            del_list.append(self.bgfmodel.a2i[i])
        del_list.sort(); del_list.reverse(); self.bgfmodel.delAtoms(del_list, False); self.bgfmodel.renumber()
        self.update()
        

    def remove_interior_water(self, *args):
        temp_list = [];
        if not self.water_position_determined:
            self.determine_water_position()

        # delete WTI atoms
        for atom in self.bgfmodel.a:
            if "I" in atom.chain and "WAT" in atom.rName:
                temp_list.append(atom.aNo)
                temp_list += atom.CONECT
        temp_list = list(set(temp_list))
        del_list = []
        for i in temp_list:
            del_list.append(self.bgfmodel.a2i[i])
        del_list.sort(); del_list.reverse(); self.bgfmodel.delAtoms(del_list, False); self.bgfmodel.renumber()
        self.update()


    def detach_hydrogen(self):
        print("Detaching hydrogens from the nanotube..")
        temp_list = []
        for atom in self.bgfmodel.a:
            if "H" in atom.ffType and "NT" in atom.rName:
                temp_list.append(atom.aNo)
        del_list = []
        for i in temp_list:
            del_list.append(self.bgfmodel.a2i[i])
        del_list.sort(); del_list.reverse(); self.bgfmodel.delAtoms(del_list, False); self.bgfmodel.renumber()
        self.update()


    def determine_water_position(self, *args):
        print("Marking water molecules whether internal or external..")
        if self.isRadiusCalc == False:
            nu.warn("Nanotube radius not calculated.")
            return 0;
        
        tempBGF = bgf.BgfFile()
        tempBGF = bgftools.make_periodic(tempBGF, self.pbc)
        self.make_pbc(self.pbc, self.bgfmodel.CRYSTX)

        # copy Nanotube atoms
        for atom in self.NTatoms + self.otheratoms:
            atom2 = copy.deepcopy(atom)
            tempBGF.addAtom(atom2)
        tempBGF.renumber()

        # mark rNames to Oxygen atoms
        n_oxygen_inside = 0; n_oxygen_outside = 0;

        for atom in self.WATatoms:
            if "H" in atom.ffType:
                continue;
        
            xyz = np.array([atom.x, atom.y, atom.z])
            r = np.sqrt(((xyz * self.axis - self.center * self.axis)**2).sum(axis=-1))
            h = (xyz * self.axismask).sum()

            if r < self.radius:
                # inside O
                atom.chain = "I"
                n_oxygen_inside += 1;
            else:
                # outside O
                atom.chain = "O"
                n_oxygen_outside += 1;

        print("Interior oxygens: " + str(n_oxygen_inside) + ", exterior oxygens: " + str(n_oxygen_outside))

        self.n_inside_water = n_oxygen_inside
        self.n_outside_water = n_oxygen_outside

        # write water molecules
        chain = ["I", "O"]
        for cname in chain:
            print("Passing " + cname)
            for atom in self.WATatoms:
                if cname in atom.chain and "O" in atom.ffType:
                    tempBGF.addAtom(atom)    #
                    l_Hatoms = atom.CONECT
                    for index, ano in enumerate(l_Hatoms):
                        atom2 = self.bgfmodel.getAtom(ano)
                        atom2.chain = cname
                        tempBGF.addAtom(atom2)    #
        self.bgfmodel = tempBGF
        self.bgfmodel.renumber()
        self.bgfmodel = bgftools.renumberMolecules(self.bgfmodel, 0, False)
        self.update()
        print("Water molecules are identified into interior(I) and exterior(O).")

        if len(args) == 0:
            nu.warn("Residue name and coords for water molecules are modified.")
        elif len(args) == 1:
            print("Modified BGF file is saved to " + str(args[0]))
            self.bgfmodel.saveBGF(args[0])

        self._water_position_determined = True


    def adjust_interior_water(self, n_water):
        
        if not self._water_position_determined:
            nu.warn("You have to run determine_water_position() first.")
            return 

        if not n_water:
            nu.warn("You must specify the number of water molecules inside the Nanotube.")
            return

        self.bgfmodel = bgftools.renumberMolecules(self.bgfmodel, 0, False)

        waters = set()
        for atom in self.bgfmodel.a:
            if atom.chain == 'I':
                waters.add(atom.rNo)
        waters = list(waters)
        print("\tFound %d interior water molecules." % len(waters))
        if len(waters) < n_water:
            nu.warn("Too few solvent molecules to choose %d molecules." % n_water)
            return

        rNos = random.sample(waters, n_water)
        print("%d water molecules are randomly chosen. Others are removed." % n_water)
        delist = []
        for atom in self.bgfmodel.a:
            if atom.chain == 'I' and not atom.rNo in rNos:
                delist.append(self.bgfmodel.a2i[atom.aNo])
        delist.sort()
        self.bgfmodel.delAtoms(delist, False)
        self.bgfmodel.renumber()
        self.bgfmodel = bgftools.renumberMolecules(self.bgfmodel, 0, False)
        self.update()
        #self.list()

        # end of adjust_interior_water()


    def list(self):
        print("\t* Height: " + str(self.height))
        print("\t* Radius: " + str(self.radius))
        print("\t* # water: " + str(len(self.WATatoms)))
        print("\t\t- inside: " + str(self.n_inside_water))
        print("\t\t- outside: " + str(self.n_outside_water))


    def add_ions_inside_nanotube(self, n_ion):
        """
        Add Na+ and Cl- ions inside nanotube.
        """
        if not self.isRadiusCalc == True:
            nu.warn("Nanotube radius not defined. Exiting.")
            return 0

        ### find the last atom before resname WAT
        aNo_lastatom = 0;
        for i in self.bgfmodel.a:
            aNo = i.aNo
            if not "WAT" in i.rName:
                if aNo > aNo_lastatom:
                    aNo_lastatom = aNo

        ### add
        n_add = 0
        while n_add < n_ion:
            x = random.uniform(0, self.pbc[0])
            y = random.uniform(0, self.pbc[1])
            z = random.uniform(0, self.pbc[2])

            xyz = np.array([x, y, z])
            r = np.sqrt(((xyz * self.axis - self.center * self.axis)**2).sum(axis=-1))
            h = (xyz * self.axismask).sum()

            if r < self.radius and self.zlo < h < self.zhi:
            
                # Na+ ion
                atom_Na = bgf.BgfAtom()
                atom_Na.x = x
                atom_Na.y = y
                atom_Na.z = z
                atom_Na.aTag = 1    # HETATM
                atom_Na.ffType = "Na"
                atom_Na.rName = "ION"
                atom_Na.aName = "Na"
                atom_Na.charge = 1.00
                atom_Na.rNo = 0
                atom_Na.chain = "X"
                
                # Cl- ion
                atom_Cl = bgf.BgfAtom()
                atom_Cl.x = x + 2
                atom_Cl.y = y + 2
                atom_Cl.z = z + 2
                atom_Cl.aTag = 1    # HETATM
                atom_Cl.ffType = "Cl"
                atom_Cl.rName = "ION"
                atom_Cl.aName = "Cl"
                atom_Cl.charge = -1.00
                atom_Cl.rNo = 0
                atom_Cl.chain = "X"

                self.bgfmodel.addAtom(atom_Na, self.bgfmodel.a2i[aNo_lastatom + 1])
                self.bgfmodel.addAtom(atom_Cl, self.bgfmodel.a2i[aNo_lastatom + 2])

                n_add += 1;

        self.bgfmodel.renumber()
        print("%d atoms added to the nanotube." % n_add)


    def evacuate_interior_water(self, *args, **kwargs):
        """
        Evacuates interior water to outside of the CNT.
        This function must be applied only to the initial structure.
        """

        # find water inside CNT
        self.determine_water_position()

        if 'delete' in kwargs and kwargs['delete'] == True:
            print("The script will remove interior water.")
            self.remove_interior_water()
        elif 'move' in kwargs and kwargs['move'] == True:
            for atom in self.WATatoms:
                if "I" in atom.chain and "O" in atom.ffType and "WAT" in atom.rName:
                    rx = self.pbc[0] - 2.0 * self.radius
                    ry = self.pbc[1] - 2.0 * self.radius
                    rz = self.pbc[2] - 2.0 * self.radius
                    r = [rx, ry, rz]

                    dx = np.random.uniform(0, rx)
                    dy = np.random.uniform(0, ry)
                    dz = np.random.uniform(0, rz)
                    d = [dx, dy, dz]
                    for index, i in enumerate(d):
                        if i > r[index] * 0.5:
                            d[index] = i + (2 * self.radius)
                    dx, dy, dz = d

                    atomH1 = self.bgfmodel.getAtom(atom.CONECT[0])
                    dx1 = atomH1.x - atom.x
                    dy1 = atomH1.y - atom.y
                    dz1 = atomH1.z - atom.z

                    atomH2 = self.bgfmodel.getAtom(atom.CONECT[1])
                    dx2 = atomH2.x - atom.x
                    dy2 = atomH2.y - atom.y
                    dz2 = atomH2.z - atom.z

                    if self.orientation == "x":
                        atom.y = dy
                        atom.z = dz
                        atomH1.y = dy + dy1
                        atomH1.z = dz + dz1
                        atomH2.y = dy + dy2
                        atomH2.z = dz + dz2

                    if self.orientation == "y":
                        atom.x = dx
                        atom.z = dz
                        atomH1.x = dx + dx1
                        atomH1.z = dz + dz1
                        atomH2.x = dx + dx2
                        atomH2.z = dz + dz2

                    if self.orientation == "z":
                        atom.x = dx
                        atom.y = dy
                        atomH1.x = dx + dx1
                        atomH1.y = dy + dy1
                        atomH2.x = dx + dx2
                        atomH2.y = dy + dy2

                    atom.chain = 'M'
                    atomH1.chain = 'M'
                    atomH2.chain = 'M'
        ### ????


    def mark_cnt_end(self):
        """
        find two ends of CNT and mark atoms' chain as E
        """

        for atom in self.NTatoms:
            if len(atom.CONECT) != 3:
                atom.chain = "E"
        print('done')

    def save(self, filename):
        if filename:
            self.bgfmodel.saveBGF(filename)
            print("Nanotube model is saved to %s " % filename)
        else:
            self.bgfmodel.saveBGF(self.model_filename)
            print("Nanotube model is saved to %s " % self.model_filename)
