# -*- coding: utf-8 -*-
"""
assign_FF.py

These classes implement a conversion utility from tinker to dl_poly.
Its primary use is to read a tinker xyz file and a tinker paramter file
and to generate the dlpoly data file from it.
It is based on the FF.py class that converts the tinker type parameter file
(key or prm) into python data structures

We use a three character string encoding for atom types. Thus, the integer numbers used in MM3
are taken literally as charcters. Any other encoding (charcters in cvff etc can also be used)

Note: Periodicity info is read from the first line after the number of atoms (as it comes from weaver!)
      Thus we currently do not take this info from the key file (it should vanish there in the pytinker framework as well)
"""
import sys
import elements
import numpy as num
import string
import copy
import itertools
from FF import FF
import numpy

from math import *
from vectortools import rotate_random
import unit_cell

# conversion factors
rad2deg = 180.0/num.pi
sigma2rmin = 2.0**(1.0/6.0)

# just needed for centering guest molecules 
atomicmass = {\
    "h" : 1.0079 , "he": 4.0026,\
    "b" : 10.811, "c" : 12.011, "n": 14.007, "o": 15.999, "f": 19.0,\
    "ne": 20.0, "na": 22.989, "s": 32, "cl": 35, "ar":40.0, "si": 28.0855, "ga":69.723, "zn": 65.409, \
    "br": 79.9,"kr": 83.0, "xe":131.0 , "cu":63.546, "se": 90.0, "zr": 91.22\
    }

# sort algorithms for types of bonds, angles, dihedrals and oops

def sort_bond(types):
    tl = string.split(types, ':')
    if cmp(tl[0], tl[1])>0:
        type = tl[1]+":"+tl[0]
    else:
        type = tl[0]+":"+tl[1]
    return type

def sort_angle(types):
    tl = string.split(types, ':')
    ta1 = tl[0]
    ta2 = tl[2]
    tc  = tl[1]
    if cmp(ta1, ta2) > 0:
        type = string.join([ta2,tc,ta1],":")
    else:
        type = string.join([ta1,tc,ta2],":")
    return type

def sort_dihedral(types):
    tl = string.split(types, ':')
    if cmp(tl[1], tl[2]) > 0:
        tl.reverse()
    if cmp(tl[1], tl[2]) == 0:
        if cmp(tl[0], tl[3]) > 0:
            tl.reverse()
    type = string.join(tl,":")
    return type

def sort_oop(types):
    tl = string.split(types, ':')
    batl = tl[1:4]
    batl.sort()
    tl[1:4] = batl[:]
    type = string.join(tl,":")
    return type


class bond:
    def __init__(self, a1, a2, types):
        if a2 > a1 :
            self.atoms = (a1, a2)
        else:
            self.atoms = (a2, a1)
        self.set_type(types[a1], types[a2])
        self.smallring=0
        self.used = False
        return
        
    def __cmp__(self, other):
        try:
            oa = other.atoms
        except AttributeError:
            a1 = other[0]
            a2 = other[1]
            if a2 > a1:
                oa = (a1, a2)
            else:
                oa = (a2, a1)
        return cmp(self.atoms, oa)
        
    def set_type(self, t1, t2):
        if cmp(t1, t2)>0:
            self.type = t2+":"+t1
        else:
            self.type = t1+":"+t2
        return
        
    
    def __repr__(self):
        return "bond %3d %3d (type %s) smallring=%d" % (self.atoms[0], self.atoms[1], self.type, self.smallring)



class angle:
    def __init__(self, aa1, ca, aa2, types):
        if aa2 > aa1:
            self.atoms =(aa1, ca, aa2)
        else:
            self.atoms =(aa2, ca, aa1)
        # set the types
        ta1 = types[aa1]
        ta2 = types[aa2]
        tc  = types[ca]
        if cmp(ta1, ta2) > 0:
            self.type = string.join([ta2,tc,ta1],":")
#            self.type_reversed = True
            if types[self.atoms[0]] != ta2:
                self.type_reversed = True
            else:
                self.type_reversed = False
        else:
            self.type = string.join([ta1,tc,ta2],":")
            self.type_reversed = False
        self.smallring = 0
        self.used = False
        return        

    def __cmp__(self, other):
        try:
            oa = other.atoms
        except AttributeError:
            aa1 = other[0]
            ca  = other[1]
            aa2 = other[2]
            if aa2 > aa1:
                oa = (aa1, ca, aa2)
            else:
                oa = (aa2, ca, aa1)
        return cmp(self.atoms, oa)
                   
    def __repr__(self):
        return "angle %d %d %d (type %s)" % (self.atoms+ (self.type,))

class dihedral:
    def __init__(self, a1, a2, a3, a4, types,smallring=0):
        if (a2 > a3):
            self.atoms = (a4, a3, a2, a1)
        else:
            if (a2 == a3):
                if (a1 > a4):
                    self.atoms = (a4, a3, a2, a1)
                else:
                    self.atoms = (a1, a2, a3, a4)
            else:
                self.atoms = (a1, a2, a3, a4)
        # set the types
        tl = map(lambda a: types[a], self.atoms)
        tl_orig = copy.deepcopy(tl)
        if cmp(tl[1], tl[2]) > 0:
            tl.reverse()
        if cmp(tl[1], tl[2]) == 0:
            if cmp(tl[0], tl[3]) > 0:
                tl.reverse()
        self.type = string.join(tl,":")
        self.type_orig = string.join(tl_orig,":")
        self.smallring = smallring
        self.used = False
        return
    
    def __repr__(self):
        return "dihedral %d %d %d %d (type %s)" % (self.atoms+ (self.type,))

    def __cmp__(self, other):
        try:
            oa = other.atoms
        except AttributeError:
            a1 = other[0]
            a2 = other[1]
            a3 = other[2]
            a4 = other[3]
            if (a2 > a3):
                oa = (a4, a3, a2, a1)
            else:
                if (a2 == a3):
                    if (a1 > a4):
                        oa = (a4, a3, a2, a1)
                    else:
                        oa = (a1, a2, a3, a4)
                else:
                    oa = (a1, a2, a3, a4)
        return cmp(self.atoms, oa)



class oop:
    # aa: central atom, ta/ba1/2 other atoms
    def __init__(self, aa, ta, ba1, ba2, types):
        oa = [ta, ba1, ba2]
        oa.sort()
        self.atoms = [aa] + oa
        tl = map(lambda a: types[a], self.atoms)
        batl = tl[1:4]
        batl.sort()
        tl[1:4] = batl[:]
        self.type = string.join(tl,":")
        self.used = False
        self.tink = False
        return
        
    def __repr__(self):
        return "oop %d %d %d %d (type %s)" % (self.atoms+ (self.type,))
        
class mol:
    
    def __init__(self, is_master):
        self.is_master = is_master
        self.verbose = is_master
        self.natoms = 0
        self.cell   = None
        self.boundarycond = 0
        self.xyz    = []
        self.cnct   = []
        self.elems  = []
        self.charges = []
        self.sigmas = None
        self.corechg= None
        self.types  = []
        self.typedata = {}
        self.bonds  = []
        self.nbonds = 0
        self.bonddata = {}
        self.nbondtypes = 0
        self.angles = []
        self.nangles = 0
        self.angledata = {}
        self.strbnddata = {}
        self.nangletypes = 0
        self.dihedrals = []
        self.ndihedrals = 0
        self.dihedraldata = {}
        self.angangdata = {}
        self.dihedraltypes = 0
        self.oops = []
        self.noops = 0
        self.oopdata = {}
        self.ooptypes = 0
        self.gaussians_used = False
        self.corechg_used   = False
        self.spbasis_used   = False
        self.var_atoms = []
        # the following data structures are for the completely reconstructed molecules/guests system
        # previously there were two seperate things: the "guests" defined via adding guests during setup
        # (but not available during a restart etc.), which could be used for the molecules_module (lambda
        # switching) in pydlpoly and the molecules detected from connectivity to be used in rigid body
        # dynamics (only used within assign_FF to write the rigid keyword to the FIELD file)
        # it is much more logic to combine and reuse both with the option of keeping the molecules
        # in the hdf5 file for a restart
        # The procedure is a follows: first the xyz file will be read and mols will be detected upon
        # connectivity. these can be named (and made rigid) via entries in the key file.
        # Then additional molecules can be added (possible also rigid). If such a system is restarted, the
        # added info during restart will be maintained.
        self.mols = []
        self.nmols= 0
        self.whichmol = []
        self.moltypes = []
        self.molnames = []
        # this datastructure is used to exclude certain terms from the potential
        # it is simply a list of flags, one for each atom. if for a given potential term
        # all atoms are excluded (True) the term will be skipped and not written to the FIELD file
        # This mechanism is used both for rigid atoms and for QM atoms in QM/MM
        self.exclude = []
        # this list contains all frozen atom indices (defined via the "freeze" keyword in the key file)
        # all frozen are excluded automatically
        self.frozen = []
        #
        self.virtual_atoms = []
        self.virtual_bonds = []
        self.virt_atom_defs = []
        self.virtual_atom_types = []
        #warning message for missing ff parameters
        self.warning_FF = False
        return
    
    def read_tinker_xyz(self, fname):
        self.xyz_filename = fname
        f = open(fname, "r")
        lbuffer = string.split(f.readline())
        self.natoms = string.atoi(lbuffer[0])
        if len(lbuffer) > 1:
            self.boundarycond = 3
            if lbuffer[1] == "#":
                # read full cellvectors
                celllist = map(string.atof,lbuffer[2:11])
                self.cell = num.array(celllist)
                self.cell.shape = (3,3)
                self.cellparams = unit_cell.abc_from_vectors(self.cell)
            else:
                try:
                    self.cellparams = map(string.atof, lbuffer[1:7])
                    self.cell = unit_cell.vectors_from_abc(self.cellparams)
                except ValueError:
                    raise ValueError('Unhashed comment after numatoms in the key file.\n...Was your key file converted from txyz format?')
            if ((self.cellparams[3]==90.0) and (self.cellparams[4]==90.0) and (self.cellparams[5]==90.0)):
                self.boundarycond=2
                if ((self.cellparams[0]==self.cellparams[1])and(self.cellparams[1]==self.cellparams[2])and\
                    (self.cellparams[0]==self.cellparams[2])):
                        self.boundarycond=1
        # intialize lists in case we reread
        self.xyz   = []
        self.elems = []
        self.types = []
        self.cnct  = []
        for i in xrange(self.natoms):
            lbuffer = string.split(f.readline())
            self.xyz.append(map(string.atof, lbuffer[2:5]))
            self.elems.append(string.lower(lbuffer[1]))
            t = lbuffer[5]
            self.types.append(t)
            if self.elems[-1] == 'xx':
                self.virtual_atoms.append(i)
                self.virtual_atom_types.append(t)
            if not self.typedata.has_key(t): self.typedata[t] = None
            self.cnct.append((num.array(map(string.atoi, lbuffer[6:]))-1).tolist())
        self.ntypes = len(self.typedata.keys())
        # check the connectivity list for virtual atoms
        self.cnct_orig = copy.deepcopy(self.cnct)
        for i in range(self.natoms):
            if i in self.virtual_atoms:
                self.virt_atom_defs.append(self.cnct[i])
                self.cnct[i] = []
            else:
                pops = []
                for j in range(len(self.cnct[i])):
                    if self.cnct[i][j] in self.virtual_atoms:
                        pops.append(self.cnct[i][j])
                       # self.cnct[i].pop(j)
                for j in pops:
                    self.cnct[i].pop(self.cnct[i].index(j))
                    #for k in pops:
                    #    self.cnct[i].pop(self.cnct[i].index(k))
        if self.verbose: print  ("$$ -- read in structure and atomtypes")
        # at this point we need to find the molecular subsystems before adding further molecules
        # --> call find_molecules at this point automatically
        self.find_molecules()
        # all atoms are included
        self.exclude = self.natoms*[False]
        # all atoms are nunfrozen
        self.frozen = self.natoms*[0]

        if "#atomtypes" in f.readline():
            stopflag = False
            self.typedata = {}
            while not stopflag:
                atmtype = f.readline().split()
                if atmtype:
                     for i in range(len(self.types)):
                         if atmtype[0] == self.types[i]:
                             if not self.typedata.has_key(atmtype[1]): self.typedata[atmtype[1]] = atmtype[2]
                             self.types[i] = atmtype[1]
                else:
                    stopflag = True
        return

    def empty_box(self, box):
        self.xyz_filename = "empty_box"
        if len(box) == 1:
            # cubic system
            self.boundarycond = 1
            self.cellparams = [box, box, box, 90.0, 90.0, 90.0]
        elif len(box) == 3:
            # orthorombic system
            self.boundarycond = 2
            self.cellparams = box + [90.0, 90.0, 90.0]
        elif len(box) == 6:
            # triclinic system
            self.cellparams = box
        else:
            raise IOError, "This type of box definition is not allowed: %s" % box
        self.cell = unit_cell.vectors_from_abc(self.cellparams)
        # set up empty data structures
        self.natoms = 0
        self.xyz   = []
        self.elems = []
        self.types = []
        self.cnct  = []
        self.ntypes = 0
        self.exclude = []
        self.frozen = []
        return
        
    def find_molecules(self):
        self.mols = []
        self.moltypes = []
        # the default moleculename is "xyz" -> molecules from the xyz file
        self.molnames = ["xyz"]
        atoms = range(self.natoms)
        self.whichmol = self.natoms * [0]
        nmol = 0
        while len(atoms) > 0:
            # localize one molecule starting from the first available atom
            leafs = [atoms[0]]
            curr_mol = []
            while len(leafs) > 0:
                new_leafs = []
                # add all to curr_mol, remove from atoms and generate new_leafs
                for l in leafs:
                    atoms.remove(l)
                    curr_mol.append(l)
                    new_leafs += self.cnct[l]
                # first remove duplicates in new_leafs
                for l in copy.copy(new_leafs):
                    i = new_leafs.count(l)
                    if i > 1:
                        for j in xrange(i-1):
                            new_leafs.remove(l)
                # now cut new_leafs (remove all those we already have in curr_mol)
                for l in copy.copy(new_leafs):
                    if curr_mol.count(l): new_leafs.remove(l)
                # now make new_leafs to leafs and continue
                leafs = new_leafs
            # at this point the molecule is complete
            curr_mol.sort()
            self.mols.append(curr_mol)
            for i in curr_mol: self.whichmol[i] = nmol
            # at this point all molecules found get the type 0 = "xyz"
            self.moltypes.append(0)
            nmol += 1
        # all atoms are assigned
        if self.verbose: print "$$ -- found %d independent molecules from connectivity" % nmol
        self.nmols=nmol
        return

    def add_mol_tinker_xyz(self, fname, N, molname, offset, scale, rotate=True):
        if self.cell is None:
            if self.verbose: print  ("ABORT ---> no periodic host system read!!")
            return
        moltype = len(self.molnames)
        self.molnames.append(molname)
        if self.verbose: print  ("$$ -- reading molecule %s from file %s" % (molname, fname))
        #if self.verbose: print  ("$$ -- WARNING this is currently limited to orthorombic systems")
        # read it in first into some temp arrays
        f = open(fname, "r")
        lbuffer = string.split(f.readline())
        mna = string.atoi(lbuffer[0])
        mxyz = []
        melems = []
        mtypes = []
        mcnct = []
        for i in xrange(mna):
            lbuffer = string.split(f.readline())
            mxyz.append(map(string.atof, lbuffer[2:5]))
            melems.append(string.lower(lbuffer[1]))
            t = lbuffer[5]
            mtypes.append(t)
            if not self.typedata.has_key(t): self.typedata[t] = None
            mcnct.append(num.array(map(string.atoi, lbuffer[6:]))-1)
        f.close()
        self.ntypes = len(self.typedata.keys())
        # center
        amass = []
        for e in melems: amass.append(atomicmass[e])
        amass = num.array(amass)
        mxyz = num.array(mxyz)
        com = sum(mxyz*amass[:,num.newaxis],0)/sum(amass)
        mxyz -= com
        # now generate N translated/rotated copies
        if self.verbose: print  ("$$ -- generating %d copies of the molecule in the box"%N)
        cell = num.array(self.cell)
        for i in xrange(N):
            act_mxyz = mxyz.copy()
            if rotate:
                # rotate molecule by random quaternion
                act_mxyz = rotate_random(act_mxyz)
            # translate
            #act_gxyz += cell*num.random.random(3)
            rand_vect = offset+(num.random.random(3)*scale)
            act_mxyz += num.sum(cell*rand_vect,axis=0)
            # register moltype (per molecule) and to which mol the atom belongs (per atom)
            self.moltypes.append(moltype)
            current_mol = self.nmols+i
            for j in xrange(mna):
                self.xyz.append(act_mxyz[j].tolist())
                self.elems.append(melems[j])
                self.types.append(mtypes[j])
                self.cnct.append((mcnct[j]+self.natoms).tolist())
                self.whichmol.append(current_mol)
            self.mols.append(range(self.natoms,self.natoms+mna))
            self.natoms += mna
        self.nmols += N
        # update exclude list
        old_natoms = len(self.exclude)
        self.exclude += (self.natoms-old_natoms)*[False]
        self.frozen  += (self.natoms-old_natoms)*[0]        
        return

    def make_extra_mol(self, atoms, name):
        """ splits atom list off from exisiting molecules and appends them as name at the end """
        mtype = len(self.molnames)
        self.molnames.append(name)
        self.moltypes.append(mtype)
        # first remove all atoms from their current molecules
        for a in atoms:
            m = self.whichmol[a]
            self.whichmol[a] = self.nmols
            self.mols[m].remove(a)
        self.mols.append(atoms)
        self.nmols += 1
        # return the index of the new molecule (always the last one)
        return self.nmols-1
        
    def find_internals(self, do_smallring_check=True, do_lin_check=False, sqp=False):
        # now find all the internal coords
        self.do_lin_check = do_lin_check
        if self.verbose: print  ("$$ -- searching for bonds")
        self.find_bonds()
        if self.verbose: print  ("$$ -- searching for angles")
        self.find_angles()
        if self.verbose: print  ("$$ -- searching for dihedrals")
        self.find_dihedrals()
        if self.verbose: print  ("$$ -- searching for oops")
        self.find_oops(sqp)
        if self.verbose: print  ("$$ -- count number of hydrogens bonded to atom")
        self.count_hydrogen()
        if do_smallring_check:
            if self.verbose: print  ("$$ -- checking for small rings")
            self.check_smallrings()
        else:
            if self.verbose: print ("$$ -- WARNING: no smallring checking done!")
        # DEBUG
        #for b in self.bonds: print b
        return
        
    def find_bonds(self):
        self.bonds=[]
        for a1 in xrange(self.natoms):
            # check all atoms this atom is bonded to
            # only add bond if the second atoms index is above (otherwise we have it already!)
            for a2 in self.cnct[a1]:
                if a2 > a1:
                    self.bonds.append(bond(a1, a2, self.types))
        self.nbonds = len(self.bonds)
        # now set up the bonddata dictionary (use tuple becasue we have pottypes)
        for b in self.bonds:
            if not self.bonddata.has_key(b.type): self.bonddata[b.type] = (None, None)
        self.nbondtypes = len(self.bonddata.keys())
        return
   
    def find_angles(self):
        self.angles=[]
        for ca in xrange(self.natoms):
            apex_atoms = self.cnct[ca]
            naa = len(apex_atoms)
            for ia in xrange(naa):
                aa1 = apex_atoms[ia]
                other_apex_atoms = apex_atoms[ia+1:]
                for aa2 in other_apex_atoms:
                    self.angles.append(angle(aa1, ca, aa2, self.types))
        self.nangles = len(self.angles)
        # now set up the angledata dictionary (use tuple becasue we use pottypes)
        for a in self.angles:
            if not self.angledata.has_key(a.type): self.angledata[a.type] = (None, None)
        self.nangletypes = len(self.angledata.keys())
        return
        
    def find_dihedrals(self):
        self.dihedrals=[]
        for a2 in xrange(self.natoms):
            # if self.verbose: print "atom %d" % a2
            for a3 in self.cnct[a2]:
                # avoid counting central bonds twice
                if a3 > a2:
                    endatom1 = list(self.cnct[a2])
                    endatom4 = list(self.cnct[a3])
                    endatom1.remove(a3)
                    endatom4.remove(a2)
                    for a1 in endatom1:
                        con1 = list(self.cnct[a1])
                        for a4 in endatom4:
                            smallring = 0
                            if a1 == a4:
                                continue 
                            if con1.count(a4):
                                smallring = 4
                            else:
                                con4 = list(self.cnct[a4])
                                for c1 in con1:
                                    if con4.count(c1):
                                        smallring = 5
                                        break
                            if self.do_lin_check == False:
                                d = dihedral(a1,a2,a3,a4, self.types, smallring=smallring)
                                #try: 
                                #    self.dihedrals.index(dihedral(a4,a2,a3,a1, self.types, smallring=smallring))
                                #except ValueError:
                                self.dihedrals.append(d)
                            else:
                                raise valueError, "Not implemented in stabel version"
        #print self.dihedrals
        self.ndihedrals = len(self.dihedrals)
        # now set up the dihedraldata dictionary (use tuple becasue of pottype)
        for d in self.dihedrals:
            if not self.dihedraldata.has_key(d.type): self.dihedraldata[d.type] = (None, None)
        self.ndihedraltypes = len(self.dihedraldata.keys())
        return

    def find_oops(self,sqp):
        self.oops=[]
        # there are a lot of ways to find oops ... 
        # we assume that only atoms with 3 partners can be an oop center
        for ta in xrange(self.natoms):
            if (len(self.cnct[ta]) == 3):
                # ah! we have an oop
                a1, a2, a3 = tuple(self.cnct[ta])
                self.oops.append(oop(ta, a1, a2, a3, self.types))
                self.oops.append(oop(ta, a2, a1, a3, self.types))
                self.oops.append(oop(ta, a3, a1, a2, self.types))
            ### not needed anymore ###
#            if ((len(self.cnct[ta]) == 3) and (self.FF.settings.has_key('opbend-next'))):
#                l = self.cnct[ta]
#                next = []
#                for i in l:
#                    if len(self.cnct[i]) == 2:
#                        if self.cnct[i][0] != ta:
#                            next.append(self.cnct[i][0])
#                        else:
#                            next.append(self.cnct[i][1])
#                if ((len(next) == 3) and (self.types[ta] in self.FF.settings['opbend-next'])):
#                    a1,a2,a3 = tuple(next)
#                    self.oops.append(oop(ta, a1, a2, a3, self.types))
#                    self.oops.append(oop(ta, a2, a1, a3, self.types))
#                    self.oops.append(oop(ta, a3, a1, a2, self.types))
        if sqp:
            raise ValueError, "Not implemented in stable version"
            #self.find_sqps()
        self.noops = len(self.oops)
        # now set up the oopdata dictionary
        for o in self.oops:
            if not self.oopdata.has_key(o.type): self.oopdata[o.type] = None
        self.nooptypes = len(self.oopdata.keys())
        return


    def count_hydrogen(self):
        self.hcount = self.natoms*[0]
        for i in xrange(self.natoms):
            hcount = 0
            for j in self.cnct[i]:
                if self.elems[j] == "h" : hcount += 1
            self.hcount[i] = hcount
        return
            
    def check_smallrings(self):
        for d in self.dihedrals:
            if d.smallring:
                self.bonds[self.bonds.index(d.atoms[0:2])].smallring = d.smallring
                self.bonds[self.bonds.index(d.atoms[1:3])].smallring = d.smallring
                self.bonds[self.bonds.index(d.atoms[2:4])].smallring = d.smallring
                self.angles[self.angles.index(d.atoms[0:3])].smallring = d.smallring
                self.angles[self.angles.index(d.atoms[1:4])].smallring = d.smallring
        return


    def load_FF(self, ff_name):
        ''' we read a key file into a FF object in order to assign the terms
        '''
        self.FF = FF(ff_name, self.verbose)
        self.FF.read_prm(ff_name)
        return
    
    def assign_FF(self, warning_FF=False):
        self.warning_FF = warning_FF
        self.assignparams_atoms()
        self.assignparams_bonds()
        self.assignparams_angles()
        self.assignparams_dihedrals()
        self.assignparams_oops()
        # virtual atoms
#        self.check_virtuals()
        return
    
    def warn_params(self, tstr, tnum, t):
        if self.warning_FF:
            trep = tstr + " " + t.replace(":", " ") + " 0.0"*tnum + "\n"
            sys.stderr.write(trep)
        else:
            raise_str = "Parameter for " + tstr + " %s not found" % t
            raise IOError, raise_str

    def assignparams_atoms(self):
        for a in self.typedata.keys():
            if self.verbose: print  "searching for typedata for atomtype %-6s" % a
            current_el = self.elems[self.types.index(a)]
            # get the atomtype outof the parameter database of FF
            aparams, equi = self.FF.get_params("atom",a,reverse=False)
            if not aparams:
                raise IOError, "Atom definition for atomtype %s is missing!!" % a
            if self.verbose: print  "element correspondence", current_el, aparams[0]
            vdwparams, equi = self.FF.get_params("vdw",a, reverse=False, verbose=False)
            m = aparams[1]
            aparams = [aparams[0]]
            if m > 0.0:
                self.typedata[a] = [aparams, vdwparams]
                self.typedata[a][0].append(m)
            else:
                self.typedata[a] = [aparams, vdwparams]
                self.typedata[a][0].append(string.atof(elements.mass[self.typedata[a][0][0]]))
            # check if this is a virtual atomm (element = xx)
            if current_el == "xx":
                self.virtual_atom_types.append(a)
        return
        
    def assignparams_bonds(self):
        equis = {}
        for b in self.bonddata.keys():
            params, pottype, equi = self.FF.get_params_and_pottype(["bond","bondq"], b)
            if params:
                self.bonddata[b] = params, pottype
                if self.verbose: print  ("$$ -- found params for bond %s (pottype: %s)" % (b, pottype))
                if equi != None:
                    equis[b] = sort_bond(equi)
                    if self.verbose: print  (" ----> bond %s equivalent to bond %s" % (b, equi))
            else:
                # first check if one of the atom types is a virtual
                virtbond = False
                for t in b.split(":"):
                    if t in self.virtual_atom_types: virtbond = True
                if not virtbond: 
                    ###### raise IOError, "Parameter for bond %s not found" % b
                    self.warn_params("bond", 2, b)
            # probe also for 5-ring params, but do not rise an error if there are none (then we use the regular ones by default!)
            params, pottype,equi = self.FF.get_params_and_pottype(["bond5"], b)
            if params:
                self.bonddata[b+"-5ring"] = params, pottype
                if self.verbose: print  ("$$ -- found params for bond %s (pottype: %s)" % (b, pottype))
                if equi != None:
                    equis[b] = sort_bond(equi)
                    if self.verbose: print  (" ----> bond %s equivalent to bond %s" % (b, equi))
            # check if there exists a virtbond entry for this bond type
            params, equi = self.FF.get_params("virtbond", b)
            if params != None:
                if self.verbose: print ("$$ -- found a virtual bond for %s" % b)
                self.virtual_bonds.append(b)
        if equis != {}:
            for i in equis.keys():
                if equis[i] not in self.bonddata.keys():
                    self.bonddata[equis[i]] = copy.deepcopy(self.bonddata[i])
                self.bonddata.pop(i)
            for b in self.bonds:
                if b.type in equis.keys():
                    b.type = equis[b.type]
        self.bond_equis = equis
        return
        
    def assignparams_angles(self):
        equis = {}
        b_equis = []
        equisc = {}
        for a in self.angledata.keys():
            # NOTE: angle, anglef and anglef-2 are mutually exclusive, that means if for a given set of atomtypes
            #       only one specific pottype can be used
            #       However, it is possible to have angle and angle5 for the same set of atomtypes because this depends on the ring type
            #        and thus we have to keep angle5 seperate 
            params, pottype, equi = self.FF.get_params_and_pottype(\
              ["angle","anglef","anglef-2","angleq"] ,a )
            if params:
                self.angledata[a] = params, pottype
                if self.verbose: print  ("$$ -- found params for angle %s (pottype: %s)" % (a, pottype))
                if equi != None:
                    s_equi = sort_angle(equi)
                    equis[a] = s_equi
                    if equi == s_equi:
                        b_equis.append(False)
                    else:
                        b_equis.append(True)
                    if self.verbose: print  (" ----> angle %s equivalent to angle %s" % (a, equi))
            else:
                # first check if one of the atom types is a virtual
                virtangle = False
                for t in a.split(":"):
                    if t in self.virtual_atom_types: virtangle = True
                if not virtangle: 
                    ######raise IOError, "Parameter for angle %s not found" % a
                    self.warn_params("angle", 2, a)
            # Now probe for angle5 .. but no error is raised if nothing is found, in this case we go back to the default
            params, pottype, equi = self.FF.get_params_and_pottype(\
              ["angle5"] ,a )
            if params:
                self.angledata[a+"-5ring"] = params, pottype
                if self.verbose: print  ("$$ -- found params for angle %s (pottype: %s)" % (a, pottype))
                if equi != None:
                    equis[a] = sort_angle(equi)
                    if self.verbose: print  (" ----> angle %s equivalent to angle %s" % (a, equi))
            # Now check for cross terms
            param, equi = self.FF.get_params("strbnd", a)
            if param:
                if self.verbose: print  ("$$ -- found also params for strbnd cross terms for %s" % a)
                self.strbnddata[a] = param
                if equi != None:
                    equisc[a] = sort_angle(equi)
                    if self.verbose: print  (" ----> strbnd %s equivalent to strbnd %s" % (a, equi))
        if equis != {}:
            for i in equis.keys():
                if equis[i] not in self.angledata.keys():
                    self.angledata[equis[i]] = copy.deepcopy(self.angledata[i])
                self.angledata.pop(i)
            for i in range(len(self.angles)):
                if self.angles[i].type in equis.keys():
                    #self.angles[i].type = equis[self.angles[i].type]
                    j = equis.keys().index(self.angles[i].type)
                    if b_equis[j] == True:
                        if self.angles[i].type_reversed == True:
                            self.angles[i].type_reversed = False
                        else:
                            self.angles[i].type_reversed = True
                    self.angles[i].type = equis[self.angles[i].type]
        if equisc != {}:
            for i in equisc.keys():
                if equisc[i] not in self.strbnddata.keys():
                    self.strbnddata[equisc[i]] = copy.deepcopy(self.strbnddata[i])
                self.strbnddata.pop(i)
#            for a in self.angles:
#                if a.type in equisc.keys():
#                    a.type = equisc[a.type]
        return

    def assignparams_dihedrals(self):
        equis = {}
        for d in self.dihedraldata.keys():
            params, pottype, equi = self.FF.get_params_and_pottype(["torsion"], d)
            if params:  
                self.dihedraldata[d] = params, pottype,
                if self.verbose: print  ("$$ -- found params for dihedral %s (pottype %s)" % (d, pottype))
                if equi != None:
                    s_equi = sort_dihedral(equi)
                    equis[d] = sort_dihedral(s_equi)
                    if self.verbose: print  (" ----> dihedral %s equivalent to dihedral %s" % (d, s_equi))
            else:
                # first check if one of the atom types is a virtual
                virtdihed = False
                for t in d.split(":"):
                    if t in self.virtual_atom_types: virtdihed = True
                if not virtdihed: 
                    ######raise IOError, "Parameter for dihedral %s not found" % d
                    self.warn_params("torsion", 4, d)
            # probe fro 5-ring params
            params, pottype, equi = self.FF.get_params_and_pottype(["torsion5"], d)
            if params:  
                self.dihedraldata[d+"-5ring"] = params, pottype,
                if self.verbose: print  ("$$ -- found params for dihedral %s (pottype %s)" % (d, pottype))
                if equi != None:
                    equis[d] = sort_dihedral(equi)
                    if self.verbose: print  (" ----> dihedral %s equivalent to dihedral %s" % (d, equi))
            # Now check for cross terms
            params, equi = self.FF.get_params("angang",d)
            if params:
                if self.verbose:print ('$$ -- found also params for angang cross terms for dihedral %s' % (d))
                self.angangdata[d] = params
        if equis != {}:
            for i in equis.keys():
                if equis[i] not in self.dihedraldata.keys():
                    self.dihedraldata[equis[i]] = copy.deepcopy(self.dihedraldata[i])
                self.dihedraldata.pop(i)
            for d in self.dihedrals:
                if d.type in equis.keys():
                    d.type = equis[d.type]
        return
        
    def assignparams_oops(self):
        equis = {}
        for o in self.oopdata.keys():
            # this is special ... we can not reverse the order of the types but we need
            # to test for any of the permutations of the final three types
            params, pottype, equi = self.FF.get_params_and_pottype(["opbend", "opbendq"], o, permute=True)
            #params = self.FF.get_params("opbend", o, permute=True)
            if params:
                self.oopdata[o] = params, pottype
                if self.verbose: print  ("$$ -- found params for oop %s" % o)
                if equi != None:
                    equis[o] = sort_oop(equi)
                    if self.verbose: print  (" ----> oop %s equivalent to oop %s" % (o, equi))
            else:
                # first check if one of the atom types is a virtual
                virtoop = False
                #for t in a.split(":"):
                for t in o.split(":"):
                    if t in self.virtual_atom_types: virtoop = True
                if not virtoop: 
                    ######raise IOError, "Parameter for oop %s not found" % o
                    self.warn_params("opbend", 2, o)
        if equis != {}:
            for i in equis.keys():
                if equis[i] not in self.oopdata.keys():
                    self.oopdata[equis[i]] = copy.deepcopy(self.oopdata[i])
                self.oopdata.pop(i)
            for o in self.oops:
                if o.type in equis.keys():
                    o.type = equis[o.type]
        return
        
    def assign_molnames(self):
        if len(self.FF.params["molname"].keys()) == 0 :
            return
        molnamedict = self.FF.params["molname"]
        molnames    = molnamedict.keys()
        for i in xrange(self.nmols):
            mol = self.mols[i]
            am  = mol[0]
            amt = self.types[am]
            named_molname = None
            for mn in molnames:
                if molnamedict[mn].count(amt):
                    if named_molname:
                        raise ValueError, "atomtype %s appears in two different molname definitions! must be unique!" % amt
                    named_molname = mn
            if named_molname:
                # yes this is potentially a named molecule, we found one of its atoms in mtype .. now check all others
                named = True
                act_moltypes = molnamedict[named_molname]
                for am in mol[1:]:
                    amt = self.types[am]
                    if not act_moltypes.count(amt): named = False
                if named:
                    # now check if named_molname is already in self.molnames
                    if not self.molnames.count(named_molname):
                        self.molnames.append(named_molname)
                    self.moltypes[i] = self.molnames.index(named_molname)
        return

    def check_virtuals(self):
        for i,e in enumerate(self.elems):
            if e == "xx":
                self.virtual_atoms.append(i)
                self.virt_atom_defs.append([])
        return
        
    def exclude_rigid(self):
        if not self.FF.settings.has_key("rigid"): return
        rigmols = self.FF.settings["rigid"]
        for i,m in enumerate(self.mols):
            if rigmols.count(self.molnames[self.moltypes[i]]):
                # this is a rigid molecule. exclude all atoms
                for j in m: self.exclude[j] = True
        return

    def exclude_frozen(self):
        if not self.FF.settings.has_key("freeze"): return
        frzmols = self.FF.settings["freeze"]
        for i,m in enumerate(self.mols):
            if frzmols.count(self.molnames[self.moltypes[i]]):
                # this is a frozen molecule. exclude and freeze all atoms
                for j in m:
                    self.exclude[j] = True
                    self.frozen[j]  = 1
        return

        
    def exclude_mol(self, molname):
        for i,m in enumerate(self.mols):
            if self.molnames[self.moltypes[i]] == molname:
                for j in m: self.exclude[j] = True
        return

    def report_molecules(self):
        if not self.verbose: return
        print ("$$ -- reporting (and consistency checking) %d molecules" % self.nmols)
        for i in xrange(self.nmols):
            atoms = self.mols[i]
            name  = self.molnames[self.moltypes[i]]
            excluded = True
            for a in atoms:
                if self.whichmol[a] != i:
                    print ("!!!!! Inconsistency in molecule assignement")
                    print ("atom %d belongs to mol %d from self.whichmol" % (a, self.whichmol[a]))
                    print ("atoms in mol: %s" % str(atoms))
                if not self.exclude[a]: excluded = False
            if excluded:
                print ("   -- molecule %3d, name %10s, number of atoms: %d (EXCLUDED)" % (i, name, len(atoms)))
            else:
                print ("   -- molecule %3d, name %10s, number of atoms: %d" % (i, name, len(atoms)))
        print ("$$ -- ")
#        print '############'
#        print bond(1,2,self.types)
#        print dihedral(1,2,3,4,self.types)
#        print angle(1,2,3,self.types)
#        print '############'
        return
        
        
#####   bookkeeping #########

    def assign_charges(self):
        charges = self.natoms*[0.0]
        if self.FF.settings.has_key("chargetype"):
            if self.FF.settings["chargetype"][0] == "none":
                self.no_charges = True
                self.charges = num.array(charges)
                return
            elif self.FF.settings["chargetype"][0] == "gaussian":
                self.gaussians_used = True
                sigmas = self.natoms*[1.0]
            elif self.FF.settings["chargetype"][0] == "gauss_core":
                self.gaussians_used = True
                self.corechg_used   = True
                sigmas = self.natoms*[1.0]
                corechg= self.natoms*[0.0]
            elif self.FF.settings["chargetype"][0] == "gauss_sp":
                self.gaussians_used = True
                self.spbasis_used = True
                sigmas = self.natoms*[1.0]
                psigmas= self.natoms*[0.0]
            else:
                raise ValueError, "Unknown charge type  %s specified" % self.FF.settings["chargetype"][0]
        # which types need modification?
        for i in xrange(self.natoms):
            t = self.types[i]
            charge_params,equi = self.FF.get_params("charge", t, verbose=False)
            if charge_params:
                charges[i] = charge_params[0]
                if self.gaussians_used: sigmas[i] = charge_params[1]
                if self.corechg_used  : corechg[i]= charge_params[2]
                if self.spbasis_used  : psigmas[i]= charge_params[2]
                # check if a modification is available
                for j in self.cnct[i]:
                    chargemod_params, equi = self.FF.get_params("chargemod", t+":"+self.types[j],reverse=False, verbose=False)
                    if chargemod_params:
                        charges[i] = chargemod_params[0]
                        if self.gaussians_used: sigmas[i] = chargemod_params[1] 
                        if self.corechg_used  : corechg[i]= charge_params[2]
                        if self.spbasis_used  : psigmas[i]= charge_params[2]
                    chargeadd_params, equi = self.FF.get_params("chargeadd", t+":"+self.types[j],reverse=False, verbose=False) ### CHECK HERE
                    if chargeadd_params:
                        charges[i]+= chargeadd_params[0]
                        if self.gaussians_used: sigmas[i]+= chargeadd_params[1] 
                        if self.corechg_used  : corechg[i]= charge_params[2]
                        if self.spbasis_used  : psigmas[i]= charge_params[2]
            # now test for "negative types" -> atom indices
            index_charge, equi = self.FF.get_params("charge", str(-i-1), verbose=False)
            if index_charge:
                charges[i] = index_charge[0]
                if self.gaussians_used: sigmas[i] = charge_params[1]
                if self.corechg_used  : corechg[i]= charge_params[2]
                if self.spbasis_used  : psigmas[i]= charge_params[2]
        self.charges = num.array(charges)
        if self.gaussians_used:
            self.sigmas  = num.array(sigmas)
        if self.corechg_used:
            self.corechg = num.array(corechg)
        if self.spbasis_used:
            self.psigmas = num.array(psigmas)
        if self.verbose:
            if self.corechg_used:
                netcore_chg = num.sum(self.corechg)
                netvale_chg = num.sum(self.charges)
                print  ("Net charge of the system: %11.5f" % (netcore_chg+netvale_chg))
                print  ("Net core charge: %10.5f    Net valence charge: %10.5f" % (netcore_chg,netvale_chg))
            else:
                print  ("Net charge of the system: %10.5f" % num.sum(self.charges))
        return

    def generate_dlpoly_typenames(self):
        self.dlp_types_dic = {}
        counter = 1
        f = open("types.dat", "w")
        for i,t in enumerate(self.types):
            if t not in self.dlp_types_dic.keys():
                self.dlp_types_dic[t]=str.upper(self.elems[i])+str(counter)
                f.write("%20s : %6s\n" % (t,self.elems[i]+str(counter)))
                counter += 1
        f.close()
        self.dlp_types = []
        for i in xrange(self.natoms):
            t = self.types[i]
            self.dlp_types.append(self.dlp_types_dic[t])
#        for i in xrange(self.natoms):
#            t = self.types[i]
#            elen = len(self.elems[i])
#            if (string.lower(t[:elen]) == self.elems[i]):
#                # ok, type starts with element code .. just keep it
#                self.dlp_types.append(t)
#            else:
#                # no .. probably just a number so we add the element
#                self.dlp_types.append(string.upper(self.elems[i])+t)
        return

#####   write FIELD file ####

    def write_FIELD(self):
#        print '############'
#        print bond(1,2,self.types)
#        print dihedral(1,2,3,4,self.types)
#        print angle(1,2,3,self.types)
#        print '############'
#        self.FF.get_variables()
        # sanity checks
#        if self.FF.settings["opbendtype"][0] != "mmff":
#            if self.verbose: print  "opbendtype needs to be MMFF to be useful with dlpoly. check your input!!"
#            return
#        if self.FF.settings["opbendpot"][0] != "dlpoly":
#            if self.verbose: print  "opbendterm needs to be DLPOLY to be useful with dlpoly. check your input!!"
#            return
#        if self.FF.settings["strbndtype"][0] != "mmff":
#            if self.verbose: print  "strbndtype needs to be MMFF to be useful with dlpoly. check your input!!"
#            return
        # do bookkeeping before output
        self.dic_strbnd = {}
        self.dic_angang = {}
        self.lreversed = [] 
        self.assign_charges()
        self.generate_dlpoly_typenames()
        f = open("FIELD", "w")
        # currently we assume that the tinker file is always one moleucle type 
        # and one moleucle of this type ... we would have to analyze the connectivity 
        # in order to distinguish
        # maybe in the future we could handle host and guest mulecules by individual
        # key files
        # we do not use neutral groups or anything like that ... third line is always empty
        f.write("%s with %s FF file\n" % (self.xyz_filename, self.FF.name))
        f.write("UNITS kcal\n")
        f.write("\n")
        f.write("MOLECULES 1\n")
        f.write("%s\n" % self.xyz_filename)
        f.write("NUMMOLS 1\n")
        # write atoms
        f.write("ATOMS %d\n" % self.natoms)
        for i in xrange(self.natoms):
            atom_param = self.typedata[self.types[i]]
            if self.gaussians_used:
                if self.corechg_used:
                    # gauss_core
                    f.write ("   %8s  %10.4f %10.4f %10.4f %10.4f 1 %1d\n" % \
                   # (self.dlp_types[i], atom_param[0][3], self.charges[i], self.sigmas[i], self.corechg[i], self.frozen[i]))
                    (self.dlp_types[i], atom_param[0][1], self.charges[i], self.sigmas[i], self.corechg[i], self.frozen[i]))
                else:
                    if self.spbasis_used:
                        # gaussians with sp-basis
                        f.write ("   %8s  %10.4f %10.4f %10.4f %10.4f 1 %1d\n" % \
                      #  (self.dlp_types[i], atom_param[0][3], self.charges[i], self.sigmas[i], self.psigmas[i], self.frozen[i]))
                        (self.dlp_types[i], atom_param[0][1], self.charges[i], self.sigmas[i], self.psigmas[i], self.frozen[i]))
                    else:
                        # gaussians
                        f.write ("   %8s  %10.4f %10.4f %10.4f 1 %1d\n" % \
                      #  (self.dlp_types[i], atom_param[0][3], self.charges[i], self.sigmas[i], self.frozen[i]))
                        (self.dlp_types[i], atom_param[0][1], self.charges[i], self.sigmas[i], self.frozen[i]))
            else:
                # regular point charges
                f.write ("   %8s  %10.4f %10.4f 1 %1d\n" % \
                (self.dlp_types[i], atom_param[0][1], self.charges[i], self.frozen[i]))
        # write bonds ... get all nonzero force constant bond terms 
        if self.verbose: print ("$$ -- %d bond terms total" % self.nbonds)
        buffer_out = ""
        count = 0
#        print self.bonddata
        for i in xrange(self.nbonds):
            b = self.bonds[i]
#            print b
            if b.type in self.virtual_bonds:
#                # first atom should be a virtual
#                if b.atoms[0] in self.virtual_atoms:
#                    vai = self.virtual_atoms.index(b.atoms[0])
#                    ra = b.atoms[1]
#                elif b.atoms[1] in self.virtual_atoms:
#                    vai = self.virtual_atoms.index(b.atoms[1])
#                    ra = b.atoms[0]
#                else:
                 raise ValueError, "somethign went wrong with the virtual atoms"
#                self.virt_atom_defs[vai].append(ra)
            else:
                bond_param, bond_pot = self.bonddata[b.type]
                if b.smallring == 5:
                    # check if there is a corresponding term -> replace params
                    btype = b.type + "-5ring"
                    if self.bonddata.has_key(btype):
                        bond_param, bond_pot = self.bonddata[btype]
                        if (('bond5' in self.FF.variables.keys()) and (b.type in self.FF.variables['bond5'].keys())): # von Johannes
                            self.FF.variables['bond5'][b.type][1].append(i)
#                if (('bond' in self.FF.variables.keys()) and (b.type in self.FF.variables['bond'].keys())): # von Johannes
#                    self.FF.variables['bond'][b.type][1].append(count)
                if ((bond_pot in self.FF.variables.keys()) and (b.type in self.FF.variables[bond_pot].keys())): # von Johannes
                    self.FF.variables[bond_pot][b.type][1].append(count)
                    self.var_atoms += b.atoms
                bond_out = self.bondterm_formatter(b.atoms, bond_param, bond_pot)
                if bond_out:
                    buffer_out += bond_out
                    b.used = True
                    count += 1
        # check if there are any distance restraints
        distrest = self.FF.params["restrain-distance"]
        for k in distrest.keys():
            params = distrest[k]
            atoms  = num.array(map(string.atoi,k.split(":")))-1
            buffer_out += self.bondterm_formatter(atoms, params, bond_pot, restraint=True)
            count += 1
        f.write("BONDS %d\n" % count)
        f.write(buffer_out)
        if self.verbose: print ("$$ -- wrote %d bond terms" % count)
        # write angles
        if self.verbose: print ("$$ -- %d angle terms total" % self.nangles)
        buffer_out = ""
        count = 0
        for i in xrange(self.nangles):
            angle = self.angles[i]
            angle_param, angle_pot = self.angledata[angle.type]
            if angle.smallring == 5:
                # check if there is a corresponding term -> replace params
                atype = angle.type + "-5ring"
                if self.angledata.has_key(atype):
                    angle_param, angle_pot = self.angledata[atype]
                    if (('angle5' in self.FF.variables.keys()) and (angle.type in self.FF.variables['angle5'].keys())): # von Johannes
                        self.FF.variables['angle5'][angle.type][1].append(count)
            if ((angle_pot in self.FF.variables.keys()) and (angle.type in self.FF.variables[angle_pot].keys())): # von Johannes
                self.FF.variables[angle_pot][angle.type][1].append(count)
                self.var_atoms += angle.atoms
            if angle_param:
                angle_out = self.angleterm_formatter(angle.atoms, angle_param, angle_pot)
                if angle_out:
                    buffer_out += angle_out
                    angle.used = True
                    count += 1
            if self.strbnddata.has_key(angle.type):
                strbnd_param = self.strbnddata[angle.type]
                b1 = bond(angle.atoms[0],angle.atoms[1],self.types)
                b2 = bond(angle.atoms[1],angle.atoms[2],self.types)
                if b1.type in self.bond_equis.keys():
                    b1_params, b1_pot = self.bonddata[self.bond_equis[b1.type]]
                    b1_type = self.bond_equis[b1.type]
                else:
                    b1_params, b1_pot = self.bonddata[b1.type]
                    b1_type = b1.type
                if b2.type in self.bond_equis.keys():
                    b2_params, b2_pot = self.bonddata[self.bond_equis[b2.type]]
                    b2_type = self.bond_equis[b2.type]
                else:
                    b2_params, b2_pot = self.bonddata[b2.type]
                    b2_type = b2.type

#                if (('bond' in self.FF.variables.keys()) and \
#                        (b1.type in self.FF.variables['bond'].keys()) and \
#                                (1 in self.FF.variables['bond'][b1.type][0].keys())):
#                    if b1.type not in self.dic_strbnd.keys():
#                        self.dic_strbnd[b1.type]=[]
#                    self.dic_strbnd[b1.type].append((count, 4))
#                if (('bond' in self.FF.variables.keys()) and \
#                        (b2.type in self.FF.variables['bond'].keys()) and \
#                                (1 in self.FF.variables['bond'][b2.type][0].keys())):
#                    if b2.type not in self.dic_strbnd.keys():
#                        self.dic_strbnd[b2.type]=[]
#                    self.dic_strbnd[b2.type].append((count, 5))
                if ((b1_pot in self.FF.variables.keys()) and \
                        (b1_type in self.FF.variables[b1_pot].keys()) and \
                                (1 in self.FF.variables[b1_pot][b1_type][0].keys())):
                    if b1_type not in self.dic_strbnd.keys():
                        self.dic_strbnd[b1_type]=[]
                    self.dic_strbnd[b1_type].append((count, 4))
                if ((b2_pot in self.FF.variables.keys()) and \
                        (b2_type in self.FF.variables[b2_pot].keys()) and \
                                (1 in self.FF.variables[b2_pot][b2_type][0].keys())):
                    if b2_type not in self.dic_strbnd.keys():
                        self.dic_strbnd[b2_type]=[]
                    self.dic_strbnd[b2_type].append((count, 5))
                if (('bond5' in self.FF.variables.keys()) and \
                        (b1.type in self.FF.variables['bond5'].keys()) and \
                                (1 in self.FF.variables['bond5'][b1.type][0].keys())):
                    if b1.type not in self.dic_strbnd.keys():
                        self.dic_strbnd[b1.type]=[]
                    self.dic_strbnd[b1.type].append((count, 4))
                if (('bond5' in self.FF.variables.keys()) and \
                        (b2.type in self.FF.variables['bond5'].keys()) and \
                                (1 in self.FF.variables['bond5'][b2.type][0].keys())):
                    if b2.type not in self.dic_strbnd.keys():
                        self.dic_strbnd[b2.type]=[]
                    self.dic_strbnd[b2.type].append((count, 5))
                # note: we need a flag for the type reversal because for an asymmetric
                #       system force constants and reference angles need to be swapped
                strbnd_out = self.strbndterm_formatter(angle.atoms, strbnd_param, \
                     angle_param, b1_params, b2_params, angle.type_reversed, angle_pot)
                if angle.type_reversed:
                    self.lreversed.append(count)
                if (('strbnd' in self.FF.variables.keys()) and (angle.type in self.FF.variables['strbnd'].keys())): # von Johannes
                    self.FF.variables['strbnd'][angle.type][1].append(count)
                if strbnd_out:
                    buffer_out += strbnd_out
                    count += 1
        # check if there are any angle restraints
        anglerest = self.FF.params["restrain-angle"]
        for k in anglerest.keys():
            params = anglerest[k]
            atoms  = num.array(map(string.atoi,k.split(":")))-1
            buffer_out += self.angleterm_formatter(atoms, params, None, restraint=True)
            count += 1
        f.write("ANGLES %d\n" % count)
        f.write(buffer_out)
        if self.verbose: print ("$$ -- wrote %d angle terms (incl. cross terms)" % count)
        # write torsions
        if self.verbose: print ("$$ -- %d dihedral terms total" % self.ndihedrals)
        buffer_out = ""
        count = 0
        for i in xrange(self.ndihedrals):
            dihed = self.dihedrals[i]
            dihed_param, dihed_pot = self.dihedraldata[dihed.type]
            if dihed.smallring == 5:
                # check if there is a corresponding term -> replace params
                dtype = dihed.type + "-5ring"
                if self.dihedraldata.has_key(dtype):
                    dihed_param, dihed_pot = self.dihedraldata[dtype]
                    if (('torsion5' in self.FF.variables.keys()) and (dihed.type in self.FF.variables['torsion5'].keys())): # von Johannes
                        self.FF.variables['torsion5'][dihed.type][1].append(count)
            if (('torsion' in self.FF.variables.keys()) and (dihed.type in self.FF.variables['torsion'].keys())): # von Johannes
                self.FF.variables['torsion'][dihed.type][1].append(count)
                self.var_atoms += dihed.atoms
            if dihed_param:
                dihed_out = self.dihedralterm_formatter(dihed.atoms, dihed_param, dihed.smallring)
                if dihed_out:
                    buffer_out += dihed_out
                    if dihed_param[:4] != 4*[0.0]:
                        dihed.used = True
                    count += 1
            if self.angangdata.has_key(dihed.type):
                angang_param = self.angangdata[dihed.type]
#                ang1=angle(dihed.atoms[0],dihed.atoms[1],dihed.atoms[2],self.types)
                ang1=dihedral(dihed.atoms[0],dihed.atoms[1],dihed.atoms[2],dihed.atoms[3],self.types)
                print ang1
                ang2=angle(dihed.atoms[1],dihed.atoms[2],dihed.atoms[3],self.types)
                ang1_params, ang1_pot = self.angledata[ang1.type]
                ang2_params, ang2_pot = self.angledata[ang2.type]
                if ((ang1_pot in self.FF.variables.keys()) and \
                        (ang1.type in self.FF.variables[ang1_pot].keys()) and \
                                (1 in self.FF.variables[ang1_pot][ang1.type][0].keys())):
                    if ang1.type not in self.dic_angang.keys():
                        self.dic_strbnd[ang1.type]=[]
                    self.dic_strbnd[ang1.type].append((count, 1))
                if ((ang2_pot in self.FF.variables.keys()) and \
                        (ang2.type in self.FF.variables[ang2_pot].keys()) and \
                                (1 in self.FF.variables[ang2_pot][ang2.type][0].keys())):
                    if ang2.type not in self.dic_angang.keys():
                        self.dic_strbnd[ang2.type]=[]
                    self.dic_strbnd[ang2.type].append((count, 2))
                angang_out = self.angangterm_formatter(dihed.atoms, angang_param, \
                        ang1_params, ang2_params)
                if (('angang' in self.FF.variables.keys()) and (dihed.type in self.FF.variables['strbnd'].keys())): # von Johannes
                    self.FF.variables['angang'][dihed.type][1].append(count)
                if angang_out:
                    buffer_out += angang_out
                    count += 1
        # check if there are any torsion restraints
        torsionrest = self.FF.params["restrain-torsion"]
        for k in torsionrest.keys():
            params = torsionrest[k]
            atoms  = num.array(map(string.atoi,k.split(":")))-1
            buffer_out += self.dihedralterm_formatter(atoms, params, None, restraint=True)
            count += 1
        f.write("DIHEDRALS %d\n" % count)
        f.write(buffer_out)
        if self.verbose: print ("$$ -- wrote %d dihedral terms" % count)
        # write opbends
        # oops are a bit special. in dlpoly only ONE term is needed for
        #   each trigonal centre since all three terms are calculated at once
        #   this implies that all terms share a common force constant
        #   NOTE: dlpoly "averages" the energy by deviding by three
        #         since tinker calculates each term seperately we simply
        #         multiply the force cosntant by three for dlpoly
        #  in order to  print  just ONE term we keep a list for the central
        #  atoms 
        if self.verbose: print ("$$ -- %d oop terms total" % self.noops)
        buffer_out = ""
        count = 0
        center_atom_list =  []
        atom_list = []
        for i in xrange(self.noops):
            oop = self.oops[i]
            otype = string.split(oop.type, ':')
            ptypes = list(itertools.permutations(otype[1:]))
            for j in range(len(ptypes)): ptypes[j] = string.join([otype[0],string.join(ptypes[j], ':')],':')
            ca = oop.atoms[0]
            oas = oop.atoms[1:]
            ptoas = list(itertools.permutations(oas))
            poas  = []
            for j in range(6): poas.append([ca]+list(ptoas[j]))
#            if not center_atom_list.count(ca):
            if ((not oop.atoms in atom_list) and (oop.tink == False)):
                # we do not have that center yet
#                center_atom_list.append(ca)
                atom_list += poas
                oop_param, oop_pot = self.oopdata[oop.type]
                oop_out = self.oopterm_formatter(oop.atoms, oop_param, oop_pot)
                if oop_out:
                    buffer_out += oop_out
                    #if (('opbend' in self.FF.variables.keys()) and (oop.type in self.FF.variables['opbend'].keys())): # von Johannes
#                    if 'opbend' in self.FF.variables.keys():
                    if oop_pot in self.FF.variables.keys():
                        for k in ptypes:
                            if k in self.FF.variables[oop_pot].keys():
                                self.FF.variables[oop_pot][k][1].append(count)
                                break
                    count += 1
                    oop.used = True
            if oop.tink == True:
                oop_param, oop_pot = self.oopdata[oop.type]
                oop_out = self.oopterm_formatter(oop.real_sort, oop_param, oop_pot, tink = True)
                if oop_out:
                    buffer_out += oop_out
                    if oop_pot in self.FF.variables.keys():
                        for k in ptypes:
                            if k in self.FF.variables[oop_pot].keys():
                                self.FF.variables[oop_pot][k][1].append(count)
                                self.var_atoms+=oop.real_sort
                                break
                    count += 1
                    oop.used = True
    #print count, oop.type, oop.atoms
        f.write("INVERSIONS %d\n" % count)
        f.write(buffer_out)
        if self.verbose: print ("$$ -- wrote %d oop terms" % count)
        # now write rigid settings (ony if there is a rigid keyword in the FF)
        if self.FF.settings.has_key("rigid"):
            rigmols = self.FF.settings["rigid"]
            buffer_out = ""
            count = 0
            for i,mt in enumerate(self.moltypes):
                if rigmols.count(self.molnames[mt]): 
                    count +=1
                    mol = num.array(self.mols[i])
                    mol += 1
                    mol = mol.tolist()
                    buffer_out += "%3d " % len(mol)
                    nlines = ((len(mol)+1)/16)+1
                    if nlines > 1:
                        pass
                    else:
                        buffer_out += (len(mol)*"%5d ")%tuple(mol)
                        buffer_out += "\n"
            f.write("RIGID %d\n" % count)
            f.write(buffer_out)
            if self.verbose: print ("$$ -- %d rigid molecules" % count)
        #
        f.write("FINISH\n")
        # write intermolecular potential here if not switched to none
        if self.FF.settings["vdwtype"][0] != "none":
            self.setup_pair_potentials()
            buffer_out = ""
            count = 0
            for pair in self.vdwdata.keys():
                vdw_out = self.vdwterm_formatter(pair, self.vdwdata[pair])
                if vdw_out:
                    buffer_out += vdw_out
                    count += 1
            f.write("VDW %d\n" % count)
            f.write(buffer_out)
            if self.verbose: print ("$$ -- wrote %d vdw pair terms" % count)
        else:
            f.write("VDW 0\n")
        f.write("CLOSE\n")
        f.close()
        return

    def all_excluded(self, atoms):
        for a in atoms:
            if not self.exclude[a]: return False
        return True

    def bondterm_formatter(self, atoms, params, pottype, restraint=False):
        if self.all_excluded(atoms): return None
        if restraint:
            k  = 2.0*params[0]*self.FF.settings["bondunit"]
            r0 = params[1]
            return "   -hrm  %5d %5d   %10.5f %10.5f\n" % (atoms[0]+1, atoms[1]+1, k, r0)
        if pottype == 'bondq':
            k  = 2.0*params[0]*self.FF.settings["bondunit"]
#            k2 = 2.0*params[2]*self.FF.settings["bondunit"]
#            k3 = 2.0*params[3]*self.FF.settings["bondunit"]
            k2 = params[2]
            k3 = params[3]
            r0 = params[1]
            return "   aqua  %5d %5d   %10.5f %10.5f %10.5f %10.5f\n" % (atoms[0]+1, atoms[1]+1, k, r0,
                k2, k3)
        # check for third param (alpha of morse) if nonzero use morse
        if params[2] != 0.0:
            # morse bond
            if self.FF.settings["bondtype"][0][:8] != "mixmorse":
                raise IOError, "you need bondtype mixmorse to allow this"
            if self.FF.settings["bondtype"][0] == "mixmorse":
                alph = params[2]
                r0   = params[1]
                k    = params[0]
                E0 = self.FF.settings["bondunit"]*k/(alph*alph)
            elif self.FF.settings["bondtype"][0] == "mixmorse_bde":
                E0 = params[2]
                r0   = params[1]
                k    = params[0]
                alph = num.sqrt(self.FF.settings["bondunit"]*k/E0)
            else:
                raise IOError, ("unknown bondtype %s " % self.FF.settings["bondtype"][0])
            return "   mors  %5d %5d   %10.5f %10.5f %10.5f\n" % \
                        (atoms[0]+1, atoms[1]+1, E0, r0, alph)
        else:
            k  = 2.0*params[0]*self.FF.settings["bondunit"]
            # if k == 0.0: return None
            r0 = params[1]
            return "   mm3b  %5d %5d   %10.5f %10.5f\n" % (atoms[0]+1, atoms[1]+1, k, r0)
    
    def angleterm_formatter(self, atoms, params, pottype, restraint=False):
        if self.all_excluded(atoms): return None
        if restraint:
            k  = 2.0*params[0]*self.FF.settings["angleunit"]*rad2deg*rad2deg
            a0 = params[1]
            # use a harmonic potential WITHOUT exclding nonboded interactions
            return "   -hrm  %5d %5d %5d   %10.5f %10.5f\n" % (atoms[0]+1, atoms[1]+1, atoms[2]+1, k, a0)
        if pottype == "angleq":
            k = 2.0*params[0]*self.FF.settings["angleunit"]*rad2deg*rad2deg
#            k2 = 2.0*params[2]*self.FF.settings["angleunit"]*rad2deg*rad2deg
#            k3 = 2.0*params[3]*self.FF.settings["angleunit"]*rad2deg*rad2deg
            k2 = params[2]
            k3 = params[3]
            a0 = params[1]
            return "   aqua  %5d %5d %5d  %10.6f %10.6f %10.6f %10.6f\n" % (atoms[0]+1, atoms[1]+1, atoms[2]+1, k, a0, k2, k3)
        if pottype == "anglef-2":
            # this is a bit trickey: the angle unit in tinker already contains the 0.5 for the
            # harmonic term which we need to remove. in addition we need to put in the multiplicity 2
            # squared to get the force constant ... this gives a factor of 2/2^2 = 1/2
            k = 0.5*params[0]*self.FF.settings["angleunit"]*rad2deg*rad2deg
            a0 = params[1]
            return "   cos2  %5d %5d %5d  %10.5f %10.5f\n" % (atoms[0]+1, atoms[1]+1, atoms[2]+1, k, a0)
        elif pottype == "anglef":
            a0 = params[1]
            fold = params[2]
            k = 0.5*params[0]*self.FF.settings["angleunit"]*rad2deg*rad2deg/fold
            flag13 = " "
            if len(params) > 3:
                # vdw and charge scaling parameters are given
                vdw13s = params[3]
                chg13s = params[4]
                if vdw13s != chg13s:
                    raise IOError, "charge and vdw scaling in an anglef term is different. This is not implemented!"
                if vdw13s == 1.0: flag13 = "-"
            return "  %1scos   %5d %5d %5d  %10.5f %10.5f %5d\n" % (flag13, atoms[0]+1, atoms[1]+1, atoms[2]+1, k, a0, fold)
        else:
            # default is harmonic (mm3 type)
            # Note: in tinker prm/key files force constants are given with respect to degrees, dlpoly
            #        works in radians. thus we have to multiply by (180/pi)**2
            k = 2.0*params[0]*self.FF.settings["angleunit"]*rad2deg*rad2deg
            # if k == 0.0 : return None
            # number of hydrogens on the central atom (this is tinker stuff and we just have to emulate
            #         this here to be consitent for the two codes)
            # Note: obviously tinker does not take the total number of hydrogens but the number of hydrogens
            #        with the two apex atoms excluded!! Thus for example the C-C-H angle on a C-CH2-C is
            #        considered as just ONE hydrogen even though CH2 has two. This is a bit awkward
            #        and i am not sure if this intended or a bug. Anyway, we do exactly the same to
            #        be sure that tinker and dlpoly get the same force constants.
            # Note2: This is only done if the second and third refangle are non-zero
            #        otherwise always the first value is used.
            if (params[2] == 0.0) and (params[3] == 0.0):
                a0 = params[1]
            else:
                hcount = self.hcount[atoms[1]]
                if self.elems[atoms[0]] == "h": hcount-=1
                if self.elems[atoms[2]] == "h": hcount-=1
                a0 =params[1+hcount]
                #if self.hcount[atoms[1]]:
                    # print ("atom %d (type %s) with angle params %s has hcount %d" % \
                    # (atoms[1]+1, self.dlp_types[atoms[1]], str(params), hcount))
            return "   mm3a  %5d %5d %5d  %10.5f %10.5f\n" % (atoms[0]+1, atoms[1]+1, atoms[2]+1, k, a0)

    def strbndterm_formatter(self, atoms, params, aparams, b1params, b2params, reversed, a_pot):
        if self.all_excluded(atoms): return None
        # Note: in tinker prm/key files force constants are given with respect to degrees, dlpoly
        #        works in radians. thus we have to multiply by (180/pi)**2
        a = 2.0*params[2]*self.FF.settings["bondunit"]
        b = params[0]*self.FF.settings["strbndunit"]*rad2deg
        c = params[1]*self.FF.settings["strbndunit"]*rad2deg
        r1= b1params[1]
        r2= b2params[1]
        # check for reversal (only flip force constants ... reference bond length are ok)
        if reversed:
            temp = b
            b = c
            c = temp
        if a_pot != 'anglef':
            t = aparams[1]
            return "   cmps  %5d %5d %5d  %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n" % \
                           (atoms[0]+1, atoms[1]+1, atoms[2]+1, \
                           a, b, c, t, r1, r2)
        else:
            t = aparams[1]
            n = aparams[2]
            return "   -cot  %5d %5d %5d  %10.5f %10.5f %10.5f %10.5f %5d %10.5f %10.5f\n" % \
                           (atoms[0]+1, atoms[1]+1, atoms[2]+1, \
                           a, b, c, t, n, r1, r2)

    def angang_formatter(self, atoms, params, a1params, a2params):
        if self.all_excluded(atoms): return None
        k = 2.0*params[0]*self.FF.settings["angleunit"]
        a1=a1params[1]
        a2=a2params[1]
        return "   cang  %5d %5d %5d %5d  %10.5f %10.5f %10.5f\n" % (atoms[0]+1, atoms[1]+1, atoms[2]+1, atoms[3+1], k, a1, a1)
    
    def dihedralterm_formatter(self, atoms, params,smallring, restraint=False):
        if self.all_excluded(atoms): return None
        if restraint:
            k  = 2.0*params[0]*self.FF.settings["angleunit"]*rad2deg*rad2deg
            t0 = params[1]
            # use a harmonic potential WITHOUT exclding nonboded interactions
            return "   harm  %5d %5d %5d %5d   %10.5f %10.5f\n" % (atoms[0]+1, atoms[1]+1, atoms[2]+1, atoms[3]+1, k, t0)
        # we can choose between cos and cos3 
        # for simplicity we put in cos3 first but if two barriers are set to zero this could be changed
        #       to reduce the effort to cos (NOTE: dl_poly uses a different form for the barrier in the two cases)
        #       in one case 0.5*V and in the other only V
        # do it quick and dirty: we should check for the phase factor! for the twofold torsion
        # it should be 180 in tinker which is equivalent to have the minimum at 0 degrees
        # in this case the co3 formula 1-cos(2phi) is correct. if the phase would be 0 then we would
        # have to make the barrier negative!
        # just do soem sanity checks
        #
        # NOTE: smallring ... usually this flag is zero.
        #  if it is 4 or 5 the dihedral term is in a four or five membered ring
        #  consequently the 1-2 or 1-3 scaling of vdw or charge interactions supersedes the usual 1-4 term
#        if not ((params[2] == 1) and (params[5] == 2) and (params[8] == 3)):
#            if self.verbose: print ("something strange with your dihedral term %s" % str(params))
        scale = self.FF.settings["torsionunit"]*2.0
#        cou14 = self.FF.settings["chg-14-scale"]
#        vdw14 = self.FF.settings["vdw-14-scale"]
        if self.FF.settings.has_key("chg-14-scale"): 
           cou14 = self.FF.settings["chg-14-scale"]
        else:
           cou14 = 1.0
        if self.FF.settings.has_key("vdw-14-scale"): 
           vdw14 = self.FF.settings["vdw-14-scale"]
        else:
           vdw14 = 1.0
        if smallring == 4:
            cou14 = 1.0
            vdw14 = 0.0
            if self.FF.settings.has_key("chg-12-scale"): cou14 = self.FF.settings["chg-12-scale"]
            if self.FF.settings.has_key("vdw-12-scale"): vdw14 = self.FF.settings["vdw-12-scale"]
        if smallring == 5:
            cou14 = 1.0
            vdw14 = 0.0
            if self.FF.settings.has_key("chg-13-scale"): cou14 = self.FF.settings["chg-13-scale"]
            if self.FF.settings.has_key("vdw-13-scale"): vdw14 = self.FF.settings["vdw-13-scale"]
        V1 = params[0]*scale
#        if (V1 != 0.0) and (params[1] != 0.0): 
#            if self.verbose: print ("phase offset for 1fold torsion not zero!")
#        V2 = params[3]*scale
        V2 = params[1]*scale
#        if (V2 != 0.0) and (params[4] != 180.0):
#            if self.verbose: print ("phase offset for 2fold torsion not 180!")
#        V3 = params[6]*scale
        V3 = params[2]*scale
#        if (V3 != 0.0) and (params[7] != 0.0):
#            if self.verbose: print ("phase offset for 3fold torsion not zero!")
#        V4 = params[9]*scale
        V4 = params[3]*scale
#        if (V4 != 0.0) and (params[10] != 180.0):
#            if self.verbose: print ("phase offset for 4fold torsion not 180!")
        if (V1==0.0) and (V2==0.0) and (V3==0.0) and (V4==0.0): 
            # no potential so we can skip if there is no 1-4 scaling
            if (cou14==1.0) and (vdw14==1.0) : return None
        if (V4 == 0.0) :
            # RS BUGFIX  ... add one 0.0 because coul4 and vdw4 must be params 5 and 6
            result = "   cos3  %5d %5d %5d %5d  %10.5f %10.5f %10.5f %10.5f  %10.5f  %10.5f\n" % \
                (atoms[0]+1, atoms[1]+1, atoms[2]+1, atoms[3]+1, V1, V2, V3, 0.0, cou14, vdw14)
        else:
            result = "   cos4  %5d %5d %5d %5d  %10.5f %10.5f %10.5f %10.5f  %10.5f %10.5f\n" % \
                (atoms[0]+1, atoms[1]+1, atoms[2]+1, atoms[3]+1, V1, V2, V3, V4, cou14, vdw14)
        return result
          
	
    def oopterm_formatter(self, atoms, params, pottype, tink = False):
        if self.all_excluded(atoms): return None
        if pottype == 'opbend':
            # tinker uses again a six order polynomial we use the default harmonic term with a phi0 of zero
            k = 2.0*params[0]*3.0*self.FF.settings["opbendunit"]*rad2deg*rad2deg
            #a0 = float(params[1])
            a0 = params[1] * (1/rad2deg)
            if k == 0.0: return None
            if tink == False:
                return "   harm  %5d %5d %5d %5d   %10.5f  %10.5f\n" % \
                       (atoms[0]+1, atoms[1]+1, atoms[2]+1, atoms[3]+1, k, a0)
                #return "   harm  %5d %5d %5d %5d   %10.5f  %10.5f\n" % \
                #       (atoms[0]+1, atoms[1]+1, atoms[2]+1, atoms[3]+1, k, 0.0)
            if tink == True:
                #return "   tink  %5d %5d %5d %5d   %10.5f  %10.5f\n" % \
                #       (atoms[0]+1, atoms[1]+1, atoms[2]+1, atoms[3]+1, k, 0.0)
                return "   tink  %5d %5d %5d %5d   %10.5f  %10.5f\n" % \
                       (atoms[0]+1, atoms[1]+1, atoms[2]+1, atoms[3]+1, k, a0)
        if pottype == 'opbendq':
            k = 2.0*params[0]*3.0*self.FF.settings["opbendunit"]*rad2deg*rad2deg
            k2 = 2.0*params[1]*3.0*self.FF.settings["opbendunit"]*rad2deg*rad2deg
            k3 = 2.0*params[2]*3.0*self.FF.settings["opbendunit"]*rad2deg*rad2deg
            if ((k == 0.0) and (k2 == 0.0) and (k3 == 0.0)): return None
            return "   quar  %5d %5d %5d %5d   %10.5f  %10.5f  %10.5f  %10.5f\n" % \
                    (atoms[0]+1, atoms[1]+1, atoms[2]+1, atoms[3]+1, k, 0.0, k2, k3)

    
    def vdwterm_formatter(self, pair, params, FIELD = True):
        t1,t2 = string.split(pair, ":")
        dlp_t1 = self.dlp_types_dic[t1]
        dlp_t2 = self.dlp_types_dic[t2]
        et1 = self.typedata[t1][0][0]
        et2 = self.typedata[t2][0][0]
        if (string.lower(t1[:len(et1)]) != et1): t1 = string.upper(et1)+t1
        if (string.lower(t2[:len(et2)]) != et2): t2 = string.upper(et2)+t2
        if (self.FF.settings["vdwtype"][0] == "exp6_damped") or \
           (self.FF.settings["vdwtype"][0] == "buckingham"):
            if self.FF.settings["raw_buck"][0] == "on":
                A = params[0]
                B = 1.0/params[1]
                C = params[2]
                D = params[3]
            else:
                A = self.FF.settings["a-expterm"]*params[1]   # abuck*eps
                #if A == 0.0: return None
                B = params[0]/self.FF.settings["b-expterm"]   # 1/(bbuck/rvdw) (dlpoly uses 1/B)
                C = params[1]*self.FF.settings["c-expterm"]*params[0]**6 # eps*cbuck*rvdw**6
                D = 0.0
            if self.FF.settings.has_key("dispscale"):
                Cscale = self.FF.settings["dispscale"]
            else:
                Cscale = 1.0
            C *= Cscale
        if self.FF.settings["vdwtype"][0] == "exp6_damped":
            Rcut = params[0]*self.FF.settings["vdwdampfact"]
            result = "   %8s %8s  ex6d  %10.2f %10.5f %10.5f %10.5f\n" % (dlp_t1, dlp_t2, A, B, C, Rcut)
            dresult = [A*418.4, B, C*418.4, Rcut, 0.0]
        elif self.FF.settings["vdwtype"][0] == "buckingham":
            if self.FF.settings["dispersion"][0] == "d3":
                A = params[1]
                B = params[0]
                C6 = params[2]
                C8 = params[3]
                Rab = params[4]
                result = "   %8s %8s  ed3   %10.2f %10.5f %10.5f %10.5f %10.5f\n" % (dlp_t1, dlp_t2, A, B, C6, C8, Rab)
                dresult = [A,B,C6,C8,Rab]
            else:
                result = "   %8s %8s  buck  %12.5f %12.5f %12.5f %12.5f\n" % (dlp_t1, dlp_t2, A, B, C, D)
                dresult = [A,B,C,D]
        elif self.FF.settings["vdwtype"][0]== "lennard-jones":
            epsilon = params[1]
            Rab     = params[0]
            result =     "   %8s %8s  lj    %12.4f %10.5f\n" % (dlp_t1, dlp_t2, epsilon, Rab)
            dresult = [epsilon, Rab,0.0,0.0,0.0]
        else:
            raise IOError, ("Sorry! vdw-type %s is not supported!" % self.FF.settings["vdwtype"][0])  
        if FIELD:
            return result
        else:
            return dresult
        
    def setup_pair_potentials(self):
        if 'vdw' in self.FF.variables.keys():
            self.var_atoms = range(self.natoms)
        # set up a matrix (upper triangle) of pairwise interactions and calculate
        # interactions via mixing parameters. then overwrite these if explicit pair
        # potentials are available
        # we work here with the values comming from the tinker key file
        # reformating and adjusting to the dl_poly format is done during output
        # by the formatter
        radrule = self.FF.settings["radiusrule"][0]
        epsrule = self.FF.settings["epsilonrule"][0]
        if self.FF.settings["radiussize"][0] == "diameter":
            radfact = 0.5
        elif self.FF.settings["radiussize"][0] == "radius":
            radfact = 1.0
        else:
            raise IOError, "The RADIUSSIZE %s in your key-file is not supported" % self.FF.settings["radiussize"][0]
        if self.FF.settings["radiustype"][0] == "r-min":
            # if it is anything like buckingham or exp6_damped all is fine, in case of lennard-jones we need to convert to sigma
            if self.FF.settings["vdwtype"][0] == "lennard-jones":
                # to convert from rmin to sigma divide by sigma2rmin
                radfact /= sigma2rmin
        elif self.FF.settings["radiustype"][0] == "sigma":
            # now complain if we do NOT use lennard-jones
            if self.FF.settings["vdwtype"][0] != "lennard-jones":
                raise IOError, "RADIUSTYPE sigma supported only for Lennard-Jones"
        else:
            raise IOError, "The RADIUSTYPE %s in your key-file is not supported" % self.FF.settings["radiustype"][0]
        self.vdwdata = {}
        types = self.typedata.keys()
        ntypes = len(types)
        for i in xrange(ntypes):
            for j in xrange(i, ntypes):
                pair = types[i]+":"+types[j]
                # check availablility of an explicit parameter
                vdwpr, equi = self.FF.get_params("vdwpr", pair, verbose=False)
                if vdwpr:
                    self.vdwdata[pair] = vdwpr
                else:
                    # get from individual params and combine
                    par_i, equi = self.FF.get_params("vdw", types[i], reverse=False, verbose=False)
                    par_j, equi = self.FF.get_params("vdw", types[j], reverse=False, verbose=False)
                    if radrule == "arithmetic":
                        rad = radfact*(par_i[0] + par_j[0])
                    elif radrule == "geometric":
                        rad = 2.0*radfact*sqrt(par_i[0]*par_j[0])
                    else:
                        if self.verbose: print  "unknown radiusrule in your tinker key file"
                    if epsrule == "arithmetic":
                        eps = 0.5 * (par_i[1]+par_j[1])
                    elif epsrule == "geometric":
                        eps = sqrt(par_i[1]*par_j[1])
                    else:
                        if self.verbose: print  "unknown epsilonrule in your tinker key file"
                    self.vdwdata[pair] = [rad, eps]
        return
    
    def write_CONFIG(self,bcond=None, vel=False):
        f = open("CONFIG","w")
        f.write("some header\n")
        velflag = 0
        if vel: velflag = 1
        if self.cell is not None:
            if bcond: self.boundarycond = bcond
            f.write("%10d%10d%10d%20.8f\n" % (velflag, self.boundarycond, self.natoms, 0.0))
            f.write("%20.12f%20.12f%20.12f\n" % tuple(self.cell[0]))
            f.write("%20.12f%20.12f%20.12f\n" % tuple(self.cell[1]))
            f.write("%20.12f%20.12f%20.12f\n" % tuple(self.cell[2]))
        else:
            f.write("%10d%10d%10d%20.8f\n" % (velflag, 0, self.natoms, 0.0))
        for i in xrange(self.natoms):
            f.write("%-8s%10d\n" % (self.dlp_types[i], i+1))
            f.write("%20.12f%20.12f%20.12f\n" % tuple(self.xyz[i]))
            if vel:
                f.write("%20.12f%20.12f%20.12f\n" % tuple(self.vel[i]))
        f.close()
        return
        
    ###################################################################    
    #################### RASPA ########################################
    ###################################################################
    
    def gen_RASPA_atomtypes(self):
        """ the problem is that RASPA does not handle individual charges per atom_param
        thus, if a chargemod is used we need another atomtype (identical to the original)
        in addition it seems RASPA does not like numbers in the atomtypes
        so we use the descriptive name in the key file and append something if a chargemod is found
        """
        self.RASPA_pseudo_atoms = {}
        self.RASPA_types =[]
        charge_params = self.FF.params["charge"]
        chrgmod_params = self.FF.params["chargemod"]
        chrgadd_params = self.FF.params["chargeadd"] ### CHECK HERE, TBI
        for i in xrange(self.natoms):
            t = self.types[i]
            RASPA_t = self.typedata[t][0][0]+t
            # check if this  atom has the original charge
            if self.charges[i] != charge_params[t][0]:
                count = 0
                for k in chrgmod_params.keys():
                    if (string.split(k,':')[0] == t) and (chrgmod_params[k][0]==self.charges[i]):
                        RASPA_t += "_%1s" % string.ascii_uppercase[count]
                    count += 1
            # check if this is already there
            if not self.RASPA_pseudo_atoms.has_key(RASPA_t):
                self.RASPA_pseudo_atoms[RASPA_t] = [t, self.charges[i]]
            self.RASPA_types.append(RASPA_t)
        return
        
    def write_RASPA_xyz(self):
        if self.verbose: print ("Writing framework structure to framework.xyz file")
        if self.verbose: print ("  --> WARNING: only orthorombic systems are supported")
        f = open("framework.xyz","w")
        f.write("%10.5f%10.5f%10.5f\n" % tuple(self.cell.diagonal()))
        f.write("%10.5f%10.5f%10.5f\n" % (90.0, 90.0, 90.0))
        f.write("%5d\n" % self.natoms)
        for i in xrange(self.natoms):
            f.write("%-10s %12.6f %12.6f %12.6f\n" % tuple([self.RASPA_types[i]]+self.xyz[i]))
        f.close()
        return
        
    def write_RASPA_cssr(self):
        if self.verbose: print ("Writing framework structure to framework.cssr file")
        if self.verbose: print ("  --> WARNING: only orthorombic systems are supported")
        # convert to fractional coordinates
        frac = num.array(self.xyz) / num.array(self.cell.diagonal())
        f = open("framework.cssr","w")
        f.write((38*" "+"%8.3f%8.3f%8.3f\n") % tuple(self.cell.diagonal()))
        f.write((21*" "+"%8.3f%8.3f%8.3f    SPGR =  1 P 1         OPT = 1\n") \
              % (90.0, 90.0, 90.0))
        f.write("%4d%4d\n" % (self.natoms,0))
        f.write("     0 %s         : %s\n" % ("RASPA_system","RASPA_system"))
        for i in xrange(self.natoms):
            f.write("%4d %-10s  %9.6f %9.6f %9.6f %4d%4d%4d%4d%4d%4d%4d%4d %7.3f\n" %\
             tuple([i+1]+[self.RASPA_types[i]]+frac[i].tolist()+8*[0]+[0.0]))
        f.close()
        return
        

    def write_RASPA_pseudo_atoms(self, extra_lines=None):
        if self.verbose: print ("Writing charges and atom types to pseudo_atoms.def")
        if extra_lines:
            if self.verbose: print ("  --> Note: added %d extra lines for guests" % len(extra_lines))
        f = open("pseudo_atoms.def","w")
        f.write("# generated by t2raspa\n")
        npseudo = len(self.RASPA_pseudo_atoms.keys())
        if extra_lines: npseudo  += len(extra_lines)
        f.write("%5d\n" % (npseudo+1))
        f.write("#type   print   as  chem   mass   charge  polarization B-factor radii  connectivity\n")
        f.write("UNIT    no      H   H      1.0      1.0     0.0        1.0      1.0     0\n")
        for t in self.RASPA_pseudo_atoms.keys():
            tinker_t = self.RASPA_pseudo_atoms[t][0]
            print tinker_t
            charge   = self.RASPA_pseudo_atoms[t][1]
            print charge
            elem     = self.typedata[tinker_t][0][0]
            print elem
            mass     = self.typedata[tinker_t][0][1]
            print mass
            f.write("%-10s%-8s%-4s%-4s%10.4f%10.4f  0.0        1.0      1.0     0\n" % \
                  (t, "yes", elem, elem, mass, charge))
        if extra_lines != None:
            for l in extra_lines: f.write(l)
        f.close()
        return
            
    ###ADDRA
    #################### RASPA v.2 ####################################

    def gen_RASPA2_atomtypes(self):
        """ the problem is that RASPA does not handle individual charges per atom_param
        thus, if a chargemod is used we need another atomtype (identical to the original)
        in addition it seems RASPA does not like numbers in the atomtypes
        so we use the descriptive name in the key file and append something if a chargemod is found
        """
        self.RASPA_pseudo_atoms = {}
        self.RASPA_reverse_pseudo_atoms = {}
        self.RASPA_types =[]
        charge_params = self.FF.params["charge"]
        chrgmod_params = self.FF.params["chargemod"]
        for i in xrange(self.natoms):
            t = self.types[i]
            RASPA_t = self.typedata[t][0][0]+t
            # check if this  atom has the original charge
            if self.charges[i] != charge_params[t][0]:
                count = 0
                for k in chrgmod_params.keys():
                    if (string.split(k,':')[0] == t) and (chrgmod_params[k][0]==self.charges[i]):
                        RASPA_t += "_%1s" % string.ascii_uppercase[count]
                    count += 1
            # check if this is already there
            ###ADDRA
            if not self.RASPA_pseudo_atoms.has_key(RASPA_t):
                self.RASPA_pseudo_atoms[RASPA_t] = [t, self.charges[i]]
            self.RASPA_types.append(RASPA_t)
            if not self.RASPA_reverse_pseudo_atoms.has_key(t):
                self.RASPA_reverse_pseudo_atoms[t] = [RASPA_t]
            elif not RASPA_t in self.RASPA_reverse_pseudo_atoms[t]:
                self.RASPA_reverse_pseudo_atoms[t].append(RASPA_t)
            ###DDARA
        return

    def write_RASPA2_pseudo_atoms(self, extra_lines=None, radii={}, connectivity={}, version=2):
        ###MOD for new RASPA formats [R.A.]
        ###TBI: radii and connectivity
        if self.verbose: print ("Writing charges and atom types to pseudo_atoms.def")
        if extra_lines:
            if self.verbose: print ("  --> Note: added %d extra lines for guests" % len(extra_lines))
        f = open("pseudo_atoms.def","w")
        f.write("# number of pseudo atoms ### generated by t2raspa v.2\n")
        npseudo = len(self.RASPA_pseudo_atoms.keys())
        if extra_lines: npseudo  += len(extra_lines)
        f.write("%5d\n" % (npseudo+1))
        f.write("#type     print   as  chem   mass     charge    polarization B-factor radii  connectivity anisotropic anisotropic-type tinker-type\n")
        f.write("UNIT      no      H   H       1.0       1.0     0.0          1.0       1.0   0            0           relative         0\n")
        for t in self.RASPA_pseudo_atoms.keys():
            tinker_t = self.RASPA_pseudo_atoms[t][0]
            charge   = self.RASPA_pseudo_atoms[t][1]
            elem     = self.typedata[tinker_t][0][0]
            mass     = self.typedata[tinker_t][0][1]
            radius   = float(radii[tinker_t][0])        ###TBI: natively in RASPA_pseudo_atoms 
            connect  = connectivity[tinker_t][1] ###TBI: natively in RASPA_pseudo_atoms
            f.write("%-10s%-8s%-4s%-4s%10.4f%10.4f  0.0          1.0     %7.3f %1s            0           absolute         %-10s\n" % \
                  (t, "yes", elem.capitalize(), elem.capitalize(), mass, charge, radius, connect, tinker_t))
        if extra_lines != None:
            for l in extra_lines: f.write(l)
        f.close()
        return

    def write_RASPA2_framework(self):
        from itertools import product as itproduct
        keydict = {
            #"atom":     (None,                  1),
            #"vdw":      (None,                  1),
            #"charge":   (None,                  1),
            #"chargemod":(None,                  1),
            "bond":     ("MM3_BOND",            2),
            "angle":    ("MM3_BEND",            2), #!!!MM3_IN_PLANE_BEND?!
            "anglef":   ("MM3_BEND",            2), #!!!MM3_IN_PLANE_BEND?!
            "torsion":  ("MM3_DIHEDRAL",        4),
            "opbend":   ("MM3_INVERSION",       4), #!!!TBI
            "strbnd":   ("MM3_BOND_BEND_CROSS", 3), #!!!TBI
            "angang":   ("MM3_BEND_BEND_CROSS", 4), #!!!TBI
            "dipole":   ("",                    2)  #!!!-2?
        }
        keyprint = {}
        for key, param_dict in self.FF.params.items():
            if key in keydict.keys():
                count = 0
                for tinker_s, param_value in param_dict.items():
                    tinker_t = tinker_s.split(":")
                    if all(t in self.RASPA_reverse_pseudo_atoms for t in tinker_t):
                        r =  [self.RASPA_reverse_pseudo_atoms[t] for t in tinker_t]
                        for k in itproduct(*r):
                            count+=1
                            kl = list(k)
                            pl = map(str,param_value)
                            ks = [ keydict[key][0] ]
                            pt = keydict[key][1]
                            keystr = "\t".join(kl+ks+pl[:pt])+"\n"
                            if key in keyprint:
                                keyprint[key][0]+=keystr
                            #    keyprint[key].append(keystr)
                            else:
                                keyprint[key] = [keystr]
                if count: keyprint[key].append(count)
        with open("framework.def","w") as f:
            keyprint["angle"][0]+=keyprint["anglef"][0]
            keyprint["angle"][1]+=keyprint["anglef"][1]
            keyord = ("coreshell","bond", "dipole","urey","angle","opbend","torsion","improper","strstr","strbnd","angang","strtor","bndtor") 
            for k in ("coreshell",        "dipole","urey",                           "improper","strstr",         "angang","strtor","bndtor"):
                keyprint[k] = ["",0]
            f.write("#CoreShells bond  BondDipoles UreyBradley bend  inv  tors improper-torsion bond/bond bond/bend bend/bend stretch/torsion bend/torsion\n")
            keycount = [keyprint[ko][1] for ko in keyord]
            f.write("          %i   %i            %i           %i   %i   %i    %i                %i         %i        %i         %i               %i            %i\n" % tuple(keycount))
            for k in keyord:
                if keyprint[k][0]:
                    f.write("#\n")
                    f.write(keyprint[k][0])

    def write_RASPA2_force_field_mixing_rules(self):
        #import pdb;pdb.set_trace()
        vdw_s="MM3_VDW_SMOOTHED5"
        vdwprint=[]
        count = 0
        for key, labs in self.RASPA_reverse_pseudo_atoms.items():
            value = self.FF.params["vdw"][key]
            for lab in labs:
                count+=1
                vdwprint.append("\t".join(["%-5s"%lab]+[vdw_s]+map(str,value)) )
        with open("force_field_mixing_rules.def","w") as f:
            f.write("# general rule for shifted vs truncated\ntruncated\n# general rule tailcorrections\nno\n# number of defined interactions")
            f.write("\n# number of defined interactions\n%i" % count)
            f.write("\n# type interaction\n")
            f.write("\n".join(vdwprint) )
            f.write("\n# general mixing rule for Lennard-Jones\nLorentz-Berthelot")
    ###DDARA

    ###################################################################    
    #################### TINKER #######################################
    ###################################################################
            
    def write_tinker_xyz(self, fname):
        if self.verbose: print  ("Writing out current molecular system as a tinker xyz file")
        f = open(fname, "w")
        f.write("%5d %10.4f %10.4f %10.4f\n" % tuple([self.natoms]+self.cell))
        for i in xrange(self.natoms):
            line = ("%3d %-3s" + 3*"%12.6f" + " %5s") % \
               tuple([i+1]+[self.elems[i]]+ self.xyz[i] + [self.types[i]])
            conn = (num.array(self.cnct[i])+1).tolist()
            if len(conn) != 0:
                line += (len(conn)*"%6d") % tuple(conn)
            f.write("%s \n" % line)
        f.close()
        return

    ###################################################################    
    #################### pdlpio /Hdf5 interface #######################
    ###################################################################
    #
    #  the purpose of these functions is to save system info into a provided pdlp file
    #  and to read also from this instead of reading tinker files
    
    def write_system_to_pdlp(self, pdlp):
        # first convert cnct info into a table
        ctab = []
        for i, ci in enumerate(self.cnct):
            for j in ci:
                if j > i : ctab.append([i,j])
        pdlp.set_system(self.elems, self.types, ctab, self.boundarycond)
        pdlp.set_molecules(self.whichmol, self.moltypes, self.molnames)
        return
        
    def get_system_from_pdlp(self, pdlp):
        """ this function reads system and molecule info from pdlp
            for a full replacement of read tinker also restart needs to be read from pdlp"""
        # system info
        self.elems, self.types, self.boundarycond, ctab = pdlp.get_system()
        self.natoms = len(self.elems)
        # recover connectivity from ctab
        self.cnct = []
        for i in xrange(self.natoms): self.cnct.append([])
        for c in ctab:
            i,j = c
            self.cnct[i].append(j)
            self.cnct[j].append(i)
        # register types
        for t in self.types:
            if not self.typedata.has_key(t): self.typedata[t] = None
        self.ntypes = len(self.typedata.keys())
        if self.verbose: print ("$$ -- read system from pdlp file")
        # molecule info - note that only a nonredundant set of info is in pdlp and the
        # rest needs to be recovered
        # NOTE: the molecule info is not detected or set but taken as found in the pdlp file
        self.whichmol, self.moltypes, self.molnames = pdlp.get_molecules()
        self.nmols = len(self.moltypes)
        self.mols = []
        for i in xrange(self.nmols) : self.mols.append([])
        for i, m in enumerate(self.whichmol): self.mols[m].append(i)
        if self.verbose: print ("$$ -- read molecule info from pdlp file")
        # set excluded atoms
        self.exclude = self.natoms * [False]
        # set frozen atoms list
        self.frozen = self.natoms * [0]
        return
        
    def read_pdlp_xyz(self, pdlp, stage, velflag, imflag):
        self.xyz_filename = pdlp.fname
        self.xyz = pdlp.read_restart(stage, "xyz").tolist()
        if velflag:
            self.vel = pdlp.read_restart(stage, "vel").tolist()
        if self.boundarycond > 0:
            self.cell = pdlp.read_restart(stage, "cell")
            self.cellparams = unit_cell.abc_from_vectors(self.cell)
        if imgflag:
            self.imgidx = pdlp.read_restart(stage, "imgidx")
        return

