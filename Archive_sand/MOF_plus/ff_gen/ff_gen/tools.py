import numpy
import string
import copy
import math
import assign_FF
import refclass
import os
import refclass

def get_multiplicity(n1, n2):
    """
    Routine to estimate m from local topology
    (stolen from QuickFF)
    """
    if   set([n1,n2])==set([5,5]): return 4
    elif set([n1,n2])==set([4,4]): return 3
    elif set([n1,n2])==set([2,4]): return 3
    elif set([n1,n2])==set([3,4]): return 3
    elif set([n1,n2])==set([3,3]): return 2
    elif set([n1,n2])==set([2,3]): return 2
    elif set([n1,n2])==set([2,2]): return 1
    else:                          return None


def get_torsion(values, m, thresshold=5):
    '''
        Get a rest value of 0.0, 360/(2*m) or None depending on the given
        equilbrium values
        (stolen from QuickFF)
    '''
    tor = [0.0, 0.0, 0.0, 0.0]
    if m == None:
        return tor
    rv = None
    per = 360/m
    for value in values:
        x = value % per
        if abs(x)<=thresshold or abs(per-x)<thresshold:
            if rv is not None and rv!=0.0:
                #tor[m-1] = 1.0
                return tor
                #return [None, None, None, None]
            elif rv is None:
                rv = 0.0
        elif abs(x-per/2.0)<thresshold:
            if rv is not None and rv!=per/2.0:
                #tor[m-1] = 1.0
                return tor
                #return [None, None, None, None]
            elif rv is None:
                rv = per/2.0
        else:
            #tor[m-1] = 1.0
            return tor
            #return [None, None, None, None]
    if rv in multidict[m]:
        tor[m-1] = 1.0
        return tor
    else:
        tor[m-1] = -1.0
        return tor


multidict = {
        1:[180.0],
        2:[0.0,180.0],
        3:[60.0,180.0,240.0],
        4:[0.0,90.0,180.0,270.0]}

def cg_filter(idx, elements, xyz):
    cg_xyz = xyz[idx]
#    for i in idx:
#        cg_xyz.append(xyz[idx,:])

    return cg_xyz, elements


class geometry:

    def __init__(self, xyz,cell=None):
        self.xyz = xyz
        self.cell = cell
        return

    def get_distance(self,atoms):
        apex_1 = numpy.array(self.xyz[atoms[0]][:])
        apex_2 = numpy.array(self.xyz[atoms[1]][:])
        if self.cell == None:
            return numpy.linalg.norm(apex_1-apex_2)

    def get_angle(self,atoms):
        apex_1 = numpy.array(self.xyz[atoms[0]][:])
        apex_2 = numpy.array(self.xyz[atoms[2]][:])
        central = numpy.array(self.xyz[atoms[1]][:])
        r1 = apex_1 - central
        r2 = apex_2 - central
        s = numpy.dot(r1,r2)/(numpy.linalg.norm(r1)*numpy.linalg.norm(r2))
        if s < -1.0: s=-1.0 
        if s >  1.0: s=1.0 
        phi = numpy.arccos(s)
        return phi * (180.0/numpy.pi)

    def get_dihedral(self, atoms):
        apex1 = numpy.array(self.xyz[atoms[0]][:])
        apex2 = numpy.array(self.xyz[atoms[3]][:])
        central1 = numpy.array(self.xyz[atoms[1]][:])
        central2 = numpy.array(self.xyz[atoms[2]][:])
        b0 = -1.0*(central1-apex1)
        b1 = central2-central1
        b2 = apex2-central2
        n1 = numpy.cross(b0,b1)
        n2 = numpy.cross(b1,b2)
        arg = -numpy.dot(n1,n2)/(numpy.linalg.norm(n1)*numpy.linalg.norm(n2))
        if abs(1.0-arg) < 10**-14:
            arg = 1.0
        elif abs(1.0+arg) < 10**-14:
            arg = -1.0
        phi = numpy.arccos(arg)
        return phi * (180.0/numpy.pi)


class atom:

    def __init__(self, element,connatoms):
        self.element = element
        connatoms.sort()
        self.connatoms = connatoms
        self.nconns = len(self.connatoms)
        self.type = ''
        return

    def __cmp__(self, other):
        if other.element != self.element:
            return cmp(self.element, other.element)
        else:
            return cmp(self.connatoms, other.connatoms)

    def get_diff(self,other):
        diff = [0,0,0]
        if __cmp__(self,other) == 0:
            return diff 
        else:
            if other.element != self.element:
                diff[0] = 1
                diff[1] = None
                diff[2] = None
                return diff
            diff[1] = other.nconns - self.nconns
            if diff[1] != 0:
                diff[2] = None
                return diff
            #diff[2] = len(set(other.connatoms)^set(self.connatoms))/float(self.nconns)
            diff[2] = len(set(other.connatoms)^set(self.connatoms))
            return diff


    def __repr__(self):
        rep =  "atom: element = %s, type = %s, bonded atoms = " % (self.element, self.type)
        rep += (self.nconns*"%s ") % tuple(self.connatoms)
        return rep

    def set_type(self,type):
        self.type = type
        return

    def get_type(self):
        return self.type

class symmetry:

    def __init__(self, xyz, masses, elems):
        self.xyz = xyz
        self.masses = masses
        self.natoms = len(masses)
        self.elems = elems
        self.symm = []
        return

    def center(self):
        molmass = numpy.sum(self.masses)
        com = numpy.sum((self.xyz*self.masses[:,numpy.newaxis]), axis = 0)
        com /= molmass
        self.com = com
        self.xyz -= self.com
        return

    def detect_symmetry(self):
        """ this is a very simple approach just checking for 
            the primary symmtry planes xy, xz, yz
        """
        other = copy.deepcopy(self)
        other.xyz = other.xyz*numpy.array([1,1,-1],"d")
        if self.mol_equal(other): self.symm.append("xy")
        other.xyz = other.xyz*numpy.array([1,-1,-1],"d")
        if self.mol_equal(other): self.symm.append("xz")
        other.xyz = other.xyz*numpy.array([-1,-1,1],"d")
        if self.mol_equal(other): self.symm.append("yz")
        return

    def mol_equal(self, other, tol=1.0e-3):
        if self.natoms != other.natoms: return False
        for i in xrange(self.natoms):
            sxyz = self.xyz[i]
            r = other.xyz-sxyz
            d = numpy.sqrt(numpy.sum(r*r, axis=1))
            closest = numpy.argsort(d)[0]
            if d[closest] > tol: return False
            if self.elems[i] != other.elems[closest]: return False
        return True

class atomtyper:

    def __init__(self, elements, xyz, cnct):
        self.elements = []
        self.avail_e  = []
        for i in elements:
            elem = string.upper(string.split(i)[0])
            if len(elem) >1:
                elem = elem[0]+string.lower(elem[1])
            self.elements.append(elem)
            if elem not in self.avail_e:
                self.avail_e.append(elem)
        self.xyz = xyz
        self.cnct = cnct
        self.natoms = numpy.shape(self.xyz)[0]
        self.atoms = []
        self.atypes = []
        self.setup_atoms()
        return

    def search_subs(self):
        from rdkit import Chem
        import IOmod
        lheavy = []
        for i in range(self.natoms):
            if self.elements[i] != 'H':
                lheavy.append(i)
        fxyz = 'tmp.xyz'
        io = IOmod.io()
        io.write_xyz(fxyz, self.elements, self.xyz)
        fmol = fxyz.split('.')[0]+'.mol'
        os.system('obabel -ixyz %s -omol > %s' % (fxyz, fmol))
        m = Chem.MolFromMolFile(fmol)
        print Chem.MolToSmiles(m)
        mnames = []
        mnatoms = []
        mindices = []
        subs = ['c1ccccc1']
        for i in subs:
            patt = Chem.MolFromSmiles(i)
            if m.HasSubstructMatch(patt):
                mnames.append(i)
                indices = list(m.GetSubstructMatches(patt))
                nindices = []
                for j in indices:
                    nindices.append([])
                    for k in j:
                        nindices[-1].append(lheavy[k])
                mindices.append(nindices)
                mnatoms.append(len(mindices[-1][0]))
        if len(mnames) > 0:
            #return  mnames, mnatoms, mindices
            sorted = numpy.argsort(numpy.array(mnatoms))
            print sorted
            mnames   = numpy.array(mnames)[sorted].tolist().reverse()
            mnatoms  = numpy.array(mnatoms)[sorted].tolist().reverse()
            mindices = numpy.array(mindices)[sorted].tolist().reverse()
            return  mnames, mnatoms, mindices



    def setup_atoms(self):
        for i in range(self.natoms):
            self.atoms.append(atom(self.elements[i], map(self.elements.__getitem__, self.cnct[i])))
        return

    def __call__(self,rules = 2): # depending on element
        """
        0 : only elements
        1 : element+coordnumber
        2 : element+coordnumber+bonded atoms
        """
        # create a rule dictionary
        self.atypes = []
        if isinstance(rules, int):
            rules = dict(zip(self.avail_e, len(self.avail_e) * [rules]))
        self.rules = rules
        # loop over all atoms
        for i,a  in enumerate(self.atoms):
            type = self.apply_rules(a)
            a.set_type(type)
            self.atypes.append(a.get_type())
        return self.atypes

    def metacall(self, short = True):
        """
        assigns atomtypes on the basis of previous assigned atomtypes,
        uses rules = 2 as basis. Can be called iteratively. 
        if short = False, it will list atomtypes of the next atoms.
        """

        types = list(set(self.atypes))
        for i in range(len(types)):
            type = types[i]
            tatoms = numpy.where(numpy.array(self.atypes)==type)[0].tolist()
            mtypes = []
            for j, t in enumerate(tatoms):
                l = []
                for k in self.cnct[t]:
                    l.append(self.atypes[k])
                mtypes.append(tuple(numpy.sort(l)))
            used = list(set(mtypes))
            nmtypes = len(used)
            if nmtypes > 1:
                ### assign new atomtypes ###
                for l in range(len(tatoms)):
                    oldtype = self.atypes[tatoms[l]]
                    if short:
                        newtype = oldtype + str(used.index(mtypes[l])+1)
                    else:
                        add = ""
                        for m in mtypes[l]:
                            add += "(" + m + ")"
                        newtype = oldtype + add
                    ### distribute new atomtypes ###
                    self.atypes[tatoms[l]] = newtype
        return


    def apply_rules(self, atom):
        rules_iml = [0,1,2]
        try:
            rule = self.rules[atom.element]
        except KeyError:
            print 'No rule found for element %s!' % atom.element
            exit()
        if rule not in rules_iml:
            print 'Rule %s not known' % rule
        type = self.apply_rule(atom, rule)
        return type

    def apply_rule(self,atom,rule):
        if rule == 0:
            type = atom.element
        elif rule == 1:
            type = atom.element+'%s' % str(atom.nconns)
        elif rule == 2:
            #type = (atom.element+'%s'+'_'+atom.nconns*"%s")\
            #        % tuple([str(atom.nconns)]+atom.connatoms)
            type = (atom.element+'%s'++atom.nconns*("%s"))\
                    % tuple([str(atom.nconns)]+atom.connatoms)
        return type


### inspired by elements.py from Arcesio
class vdwp:

    def __init__(self, path = None):
        import csv
        import os
        from collections import OrderedDict
        params_in = []
        types_in = []
        path = os.path.join(os.environ['FFDIR'],'ff_gen/ext_params.csv')
        with open(path, 'r',) as fparam:
            reader = csv.reader(fparam, delimiter=',', skipinitialspace = True)
            reader.next()
            for row in reader:
                param = []
                types_in.append(row[0])
                param.append(float(row[1]))
                param.append(float(row[2]))
                params_in.append(param)

        self.params = OrderedDict(zip(types_in, params_in))
        return

    def __call__(self,types, i = 0):
        self.set= []
        if i > len(types)-1:
            print 'ERROR: No vdw parameter set found for element %s' % types[2]
            exit()
        try:
            self.set = self.params[types[i]]
        except KeyError:
            self.__call__(types, i+1)
            
class keycreator:

    def __init__(self, fxyz=None, fref = None):
        if fxyz != None:
            self.pass_coordinates(fxyz,fref)
#        self.fxyz = fxyz
#        self.mol = assign_FF.mol(is_master=True)
#        self.mol.verbose = False
#        self.mol.read_tinker_xyz(self.fxyz)
#        self.mol.find_internals(do_smallring_check=False, do_lin_check=False, sqp = False)
#        self.geom = tools.geometry(self.mol.xyz)
#        if fref != None:
#            refdata = [self.mol.xyz, self.mol.elems]
#            self.ref = refclass.refclass(fref, refdata)
        self.formatted_keys = {
            "bond"     : 2,
            "atom"     : 1,
            "opbend"   : 4,
            "angle"    : 3,
            "anglef"   : 3,
            "strbnd"   : 3,
            "vdw"      : 1,
            "charge"   : 1,
            "strbnd"   : 3,
            "torsion"  : 4}
        self.var_settings = {
                'opbend' : [[[0.0, 0.5],'a']],
                'torsion': [[[0.0, 5.0],'a']],
                'angle'  : [[[0.5,2.0],['a']],[[0.94,1.06],'r']],
                'strbnd' : [[[0.0,0.5],['a']]],
                'bond'   : [[[0.5,3.0],['a']],[[0.94,1.06],'r']]}
        self.head = "version             2.0\n\
parameters          none\n\
bondunit            71.94\n\
angleunit           0.02191418\n\
opbendunit          0.02191418\n\
strbndunit          2.51118\n\
torsionunit         0.5\n\
vdwtype             exp6_damped\n\
vdwdampfact         0.25\n\
radiusrule          arithmetic\n\
radiustype          r-min\n\
radiussize          radius\n\
epsilonrule         geometric\n\
a-expterm           184000.0\n\
b-expterm           12.0\n\
c-expterm           2.25\n\
bondtype            mixmorse_bde\n\
strbndtype          mmff\n\
opbendtype          mmff\n\
chargetype          gaussian\n\n"

    def pass_mol(self,mol):
        self.mol = mol
        return


    def pass_coordinates(self,fxyz,fref):
        self.fxyz = fxyz
        self.mol = assign_FF.mol(is_master=True)
        self.mol.verbose = False
        self.mol.read_tinker_xyz(self.fxyz)
        self.mol.find_internals(do_smallring_check=False, do_lin_check=False, sqp = False, 
                ltorsion = False)
        self.geom = geometry(self.mol.xyz)
        if fref != None:
            refdata = [self.mol.xyz, self.mol.elems]
            self.ref = refclass.refclass(fref, refdata)
        return

    def get_atoms(self):
        for i in self.mol.typedata.keys():
            for j,t in enumerate(self.mol.types):
                if i == t:
                    self.mol.typedata[i] = [[self.mol.elems[j]],[[],[],[],[]],[[]]]
                    break


    def get_charges(self, tag = 'primary'):
        if hasattr(self,'ref') == False:
            return
        qmcharges = numpy.around(copy.deepcopy(self.ref(info = 'charges', tag = tag).T),6)
        qmnetcharge = numpy.sum(qmcharges)
        #print qmnetcharge
        if abs(qmnetcharge) > 0.005:
            print 'WARNING: QM system seems not to be zero charged!!!'
        else:
            qmcharges[0] -= qmnetcharge
            qmnetcharge = numpy.sum(qmcharges)
        radii = {
                'h' : 0.723638,
                'c' : 1.162986,
                'zn': 2.073300,
                'n' : 1.125000,
		'ni': 2.073300,
                'o' : 1.117553,
                'na': 4.146532, # roberto (invented)
                'cu': 2.073266}
        for i,t in enumerate(self.mol.types):
            self.mol.typedata[t][1][0].append(qmcharges[i,0])
        charges = []
        for i in self.mol.typedata:
            self.mol.typedata[i][1][1]=numpy.round(numpy.mean(self.mol.typedata[i][1][0]),6)
            self.mol.typedata[i][1][2]=numpy.std(self.mol.typedata[i][1][0])
            self.mol.typedata[i][1][3]=radii[self.mol.typedata[i][0][0]]
            for j in range(len(self.mol.typedata[i][1][0])):
                charges.append(self.mol.typedata[i][1][1])
        mmnetcharge = numpy.sum(charges)
        return


    def get_vdws(self):
        type = atomtyper(self.mol.elems, self.mol.xyz, self.mol.cnct)
        t1 = copy.deepcopy(type(rules = 1))
        t2 = copy.deepcopy(type(rules = 2))
        #print t2
        #print t1
        vdw = vdwp('/home/johannes/work/FFs/CuFormate/test/test.csv')
#        for i in range(len(t1)):
#            vdw([t2[i], t1[i], self.mol.elems[i]])
#            params = vdw.set
#            print params
        for i, t in enumerate(self.mol.types):
            self.mol.typedata[t][2][0].append([t2[i],t1[i],self.mol.elems[i]])
        for i in self.mol.typedata:
            t1s = []
            t2s = []
            for j in range(len(self.mol.typedata[i][2][0])):
                t1s.append(self.mol.typedata[i][2][0][j][1])
                t2s.append(self.mol.typedata[i][2][0][j][0])
                if t2s.count(t2s[0])==len(t2s):
                    types = self.mol.typedata[i][2][0][0]
                elif t1s.count(t1s[0])==len(t1s):
                    types = self.mol.typedata[i][2][0][0][1:]
                else:
                    types = self.mol.typedata[i][2][0][0][2:]
            vdw(types)
            self.mol.typedata[i][2].append(vdw.set)

           # print self.mol.typedata[t]
        return


    def get_angles(self):
        for a in self.mol.angledata:
            self.mol.angledata[a] = [[],[[],[],[]]]
        for a in self.mol.angles:
            phi = self.geom.get_angle(a.atoms)
            self.mol.angledata[a.type][1][0].append(phi)
        for a in self.mol.angledata:
            self.mol.angledata[a][1][1]=numpy.mean(self.mol.angledata[a][1][0])
            self.mol.angledata[a][1][2]=numpy.std(self.mol.angledata[a][1][0])
        return


    def get_bonds(self):
        for a in self.mol.bonddata:
            self.mol.bonddata[a] = [[],[[],[],[]]]
        for a in self.mol.bonds:
            r= self.geom.get_distance(a.atoms)
            self.mol.bonddata[a.type][1][0].append(r)
        for a in self.mol.bonddata:
            self.mol.bonddata[a][1][1]=numpy.mean(self.mol.bonddata[a][1][0])
            self.mol.bonddata[a][1][2]=numpy.std(self.mol.bonddata[a][1][0])
        return

    def get_dihedrals(self):
        for a in self.mol.dihedraldata:
            self.mol.dihedraldata[a] = [[[],[]],[[],[]]]
        for a in self.mol.dihedrals:
            phi = self.geom.get_dihedral(a.atoms)
            m   = get_multiplicity(len(self.mol.cnct[a.atoms[1]]),len(self.mol.cnct[a.atoms[2]]))
            self.mol.dihedraldata[a.type][1][0].append(phi)
            self.mol.dihedraldata[a.type][0][0].append(m)
        for a in self.mol.dihedraldata.keys():
            if None in self.mol.dihedraldata[a][0][0]:
                self.mol.dihedraldata[a][0][1] = None
            else:
                self.mol.dihedraldata[a][0][1]=int(numpy.round(numpy.mean(self.mol.dihedraldata[a][0][0])))
            self.mol.dihedraldata[a][1][1]=get_torsion(self.mol.dihedraldata[a][1][0],
                    self.mol.dihedraldata[a][0][1],20)
            #print a 
            #print self.mol.dihedraldata[a]
            #print '--------------------------'
        return


    def write_empty_key(self,fkey, var = False, strbnd = True):
        buffer_types     = ""
        buffer_bonds     = ""
        buffer_angles    = ""
        buffer_strbnd    = ""
        buffer_dihedrals = ""
        buffer_opbends   = ""
        # nonbonded
        for i, t in enumerate(self.mol.typedata.keys()):
            buffer_types += self.formatter("atom", [t], params = self.mol.typedata[t][0])
            buffer_types += self.formatter("vdw", [t], params = self.mol.typedata[t][2][1])
            if hasattr(self,'ref'):
                buffer_types += self.formatter("charge", [t], params = [self.mol.typedata[t][1][1]]+\
                        [self.mol.typedata[t][1][3]], fieldwidth=18)
            else:
                buffer_types += self.formatter("charge", [t])
        # bonded
        for i, b in enumerate(self.mol.bonddata.keys()):
            buffer_bonds += self.formatter("bond", string.split(b,':'),var, 
                    params = [2.0]+[self.mol.bonddata[b][1][1]])
        for i, a in enumerate(self.mol.angledata.keys()):
            buffer_angles += self.formatter("angle", string.split(a,':'),var,
                    params = [1.0]+[self.mol.angledata[a][1][1]])
        if strbnd == True:
            for i, a in enumerate(self.mol.angledata.keys()):
                buffer_strbnd += self.formatter("strbnd", string.split(a,':'),var,
                        params = [0.2,0.2,0.2])
        for i, d in enumerate(self.mol.dihedraldata.keys()):
            buffer_dihedrals += self.formatter("torsion", string.split(d,':'),var, 
                    params = self.mol.dihedraldata[d][1][1])
        for i, o in enumerate(self.mol.oopdata.keys()):
            buffer_opbends += self.formatter("opbend", string.split(o,':'),var,
                    params = [0.1])
        buffer_total = self.head + buffer_types +"\n" + buffer_bonds +"\n" + buffer_angles +"\n"\
                + buffer_strbnd +"\n" + buffer_dihedrals + "\n" + buffer_opbends
        f = open(fkey, 'w')
        f.write(buffer_total)
        f.close
        return

    def buffer_nonbonded(self,):
        buffer_types = ''
        for i, t in enumerate(self.mol.typedata.keys()):
            buffer_types += self.formatter("atom", [t], params = self.mol.typedata[t][0])
            buffer_types += self.formatter("vdw", [t], params = [1.0,0.0])
            if hasattr(self,'ref'):
                buffer_types += self.formatter("charge", [t], params = [self.mol.typedata[t][1][1]]+\
                        [self.mol.typedata[t][1][3]], fieldwidth=18)
            else:
                buffer_types += self.formatter("charge", [t], params = [0.0,1.0])
        return buffer_types

    def buffer(self,keyword,attr,var = False):
        buffer = ''
        data = getattr(self.mol,attr)
        for i, b in enumerate(data.keys()):
            if ((keyword == 'angle') and (len(data[b])>2)): keyword = 'anglef'
            buffer += self.formatter(keyword, string.split(b,':'),var, 
                    params = data[b])
            if keyword == 'anglef': keyword = 'angle'
        return buffer

    def write_key(self,fkey,buffers):
        f = open(fkey, 'w')
        f.write(self.head)
        for i in buffers: f.write(i+"\n")
        f.close
        return

    def formatter(self,keyword, atypes,var=False,fieldwidth = 12,params = None):
        na = self.formatted_keys[keyword]
        form = "%-15s" + na*"%-5s"
        if  params != None:
            for i in params:
                if isinstance(i, basestring):
                    form += "%-5s"
                elif isinstance(i, float):
                    form += "%%%d.%df" % (fieldwidth, fieldwidth/3)
                elif isinstance(i, int):
                    form += "%%%dd" % fieldwidth
            formo = form % tuple([keyword]+atypes+params)
        else:
            formo = form % tuple([keyword]+atypes)
        formo += "\n"
        if var == True:
            varformo = ''
            for i in range(len(params)):
                pstrin = "%%%d.%df" % (fieldwidth, fieldwidth/3)
                varform = "var %-11s" + na*"%-5s" + "%-5s"+2*pstrin 
                varform += "\n"
                try:
                    bounds = copy.deepcopy(self.var_settings[keyword][i])
                except IndexError:
                    bounds = copy.deepcopy(self.var_settings[keyword][0])
                if bounds[1] == 'r':
                    bounds[0][0] *=params[i]
                    bounds[0][1] *=params[i]
                if keyword == 'torsion' and params[i] == 0.0:
                    pass
                else:
                   if keyword == 'torsion': ### handle negative torsions
                       shift = bounds[0][1]*(1-int(params[i]))/2 ###0 if 1, 5 if -1
                       bounds[0][0] -= shift
                       bounds[0][1] -= shift
                   varformo += varform % tuple([keyword]+atypes+[i+1]+
                        bounds[0])
            if keyword == 'strbnd':
                if atypes[0] == atypes[2]:
                    varformo += ("var con  %-11s" + na * "%-10s"  + " 1  " \
                            + na* "%-10s" + "  2  \n") % tuple([keyword] + atypes +atypes)
            return formo + varformo
        return formo


