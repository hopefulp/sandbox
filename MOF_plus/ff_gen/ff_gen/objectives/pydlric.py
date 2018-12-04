"""
This file implements a ric_fit class.
It is inheriting all features from RedIntCoords in ric.py
and adds features to load a respective reference structure and Hessian.
In addition, a weight matrix is held to evaluate various kinds of weighted mean
square deviations to be used as ingredients to fittness values
"""

import string
import numpy as np
from ff_gen.ric import RedIntCoords 
import  ff_gen.tools as tools
import copy
import os


def revlist(l):
    nl = copy.copy(l)
    nl.revert()
    return nl

ricdic = {'str':[2,0],
        'ibe' : [3,1],
        'obe' : [4,2],
        'lbe' : [3,3],
        'tor' : [12,4]}


class pydlric():
    """
    class to compute redundant internal coordinates (ric)
    by using the inherited ric module and to compute deviations from a corresponding
    reference.
    """
    
    def __init__(self):
        return
    

    def prepare(self):
        """
        Setup the RICs 
        
        Call after adding all rics
        
        :Parameters:
        
            - masses : numpy array [natoms] defining the atomic masses
            
        """
        # now let us do all the extra stuff to map rics and allocate python arrays for the reference
        self.nstretch     = len(self.pd.ric._stretches)
        self.first_str    = 0
        self.nin_bend     = len(self.pd.ric._in_bends)
        self.first_ibe    = self.first_str+self.nstretch
        self.nout_bend    = len(self.pd.ric._out_bends)
        self.first_obe    = self.first_ibe+self.nin_bend
        self.nlin_bend    = len(self.pd.ric._lin_bends)
        self.first_lbe    = self.first_obe+self.nout_bend
        self.ntorsion     = len(self.pd.ric._torsions)
        self.first_tor    = self.first_lbe+self.nlin_bend
        # keep a list of central atoms for all torsions
        self.tor_cent = []
        for t in self.pd.ric._torsions:
            na = len(t)
            self.tor_cent.append(t[na/2-1:na/2+1])
        self.cycle = 0
        self.riclist = []
        self.riclist.append(self.pd.ric._stretches)
        self.riclist.append(self.pd.ric._in_bends)
        self.riclist.append(self.pd.ric._out_bends)
        self.riclist.append(self.pd.ric._lin_bends)
        self.riclist.append(self.pd.ric._torsions)
        self.firstrics = []
        self.firstrics.append(self.first_str)
        self.firstrics.append(self.first_ibe)
        self.firstrics.append(self.first_obe)
        self.firstrics.append(self.first_lbe)
        self.firstrics.append(self.first_tor)
        return
    
    def map_ric(self, ric_type, ind):
        """
        Map a ric definition to its internal index
        
        :Parameters:
        
            - ric_type : string defining the type of ric ("str", "ibe", "obe", "lbe", "tor")
            - ind      : index of the defining atoms
            
        :Returns:
        
            - index within ric_type to be used with the value arrays
            - global index within ric Hessian
            
        :Note: 
            Since we do not know the ordering of atoms for torsions we identify torsions by the central
            two atoms. In other words there can be only one torsion ric for a rotation around a given 
            bond. 
        """
#
        assert len(ind) == ricdic[ric_type][0]
        try:
            i = self.riclist[ricdic[ric_type][1]].index(ind)
        except ValueError:
            ind.reverse()
            try:
                i = self.riclist[ricdic[ric_type][1]].index(ind)
            except ValueError:
                return -1, -1
        iglob = self.firstrics[ricdic[ric_type][1]]+i
        return i, iglob

    def map_fractor_ric(self, ind):
       assert len(ind) == 2
       try:
           i = self.tor_cent.index(ind)
       except ValueError:
           ind.reverse()
           try:
               i = self.tor_cent.index(ind)
           except ValueError:
               return -1, -1
       iglob = self.first_tor+i
       return i, iglob

    def read_lin_ref(self, fname):
        f = open(fname, 'r')
        self.lin_dict = {}
        line = f.readline()
        stop = False
        while not stop:
            sline = string.split(line)
            self.lin_dict[tuple(map(string.atoi,sline[0:3]))] = string.atoi(sline[3])
#            self.lin_dict[string.atoi(sline[0])] = string.atoi(sline[1])
            line = f.readline()
            if len(line) == 0:
                f.close
                stop = True
        return
  
    def get_ric(self, lin_autoassign=True):
        """
        Gets the RICs of the actual system
        """
        ### bonds ###
        lbonds = []
        for b in self.pd.mol.bonds:
            if b.used == True:
                lbonds.append(b.atoms)
        self.abonds = np.array(lbonds) + 1
        ### angles ###
        langles = []
        llinangles = []
        self.lin_ref_atoms = []
        lin_flag=False
        for a in self.pd.mol.angles:
            if a.used == True:
                if self.check_lin(a.atoms) == False:
                    langles.append(a.atoms)
                    a.lin = False
                if self.check_lin(a.atoms) == True:
                    print a.atoms
                    a.lin = True
                    llinangles.append(a.atoms)
                    if hasattr(self,'lin_dict'): ###read_lin_ref explicitly called
                        if a.atoms in self.lin_dict.keys():
                            self.lin_ref_atoms.append(self.lin_dict[a.atoms])
                    else:
                        if lin_autoassign:
                            self.lin_ref_atoms.append(self.choose_lin(a.atoms))
                        else:
                            lin_flag=True
                            print 'no ref atom definition found for lin bend %d' % (a.atoms[1]+1)
                            #raise IOError
        if lin_flag: raise IOError
        self.lin_dict = dict(zip(llinangles, self.lin_ref_atoms))
        self.aangles = np.array(langles) + 1
        self.alinangles = np.array(llinangles) + 1
        ### oop ###
        loops = []
        i = 0
        for o in self.pd.mol.oops:
           if o.tink == False:
                if (i/3.0 - i/3) == 0.0:
                    a = o.atoms[0]
                    b = o.atoms[1]
                    c = o.atoms[2]
                    d = o.atoms[3]
                    loops.append([a,b,c,d])
                    loops.append([a,c,b,d])
                    loops.append([a,d,b,c])
                i += 1
           else:
                pass
                #loops.append(o.real_sort)
        self.aoops = np.array(loops) +1
        ### dihedrals ###
        if self.fragtor == False:
            self.get_ric_torsions()
            return
        ldihedrals = []
        for i in range(len(lbonds)):
            a1 = lbonds[i][0]
            a2 = lbonds[i][1]
            types = []
            types.append(self.pd.get_atomtypes()[a1])
            types.append(self.pd.get_atomtypes()[a2])
            if types in self.notwgh:
                continue
            if ((len(self.pd.mol.cnct[a1]) > 1 ) and (len(self.pd.mol.cnct[a2]) > 1)): 
                ca1 = copy.deepcopy(self.pd.mol.cnct[a1])
                ca1.remove(a2)
                ca2 = copy.deepcopy(self.pd.mol.cnct[a2])
                ca2.remove(a1)
                # check for linear/undefined dihedral angles
                if len(ca1) == 1:
                    if self.check_lin([a2,a1,ca1[0]]) == True:
                        continue
                if len(ca2) == 1:
                    if self.check_lin([a1,a2,ca2[0]]) == True:
                        continue
                a1 += 1
                a2 += 1
                ca1 = list(np.array(ca1)+1)
                ca2 = list(np.array(ca2)+1)
                for j in ca1:
                    if j in ca2:
                        ca2.pop(ca2.index(j))
                if len(ca1) == len(ca2):
                    ldihedrals.append(ca1+[a1]+[a2]+ca2)
                else:
                    if len(ca1) > len(ca2):
                        nzeros = len(ca1) - len(ca2)
                        ca2 += nzeros * [0]
                    elif len(ca1) < len(ca2):
                        nzeros  = len(ca2) - len(ca1)
                        ca1 += nzeros * [0]
                    ldihedrals.append(ca1+[a1]+[a2]+ca2)
        self.adihedrals = ldihedrals
        return

    def get_ric_torsions(self):
        ldihedrals = []
        for d in self.pd.mol.dihedrals:
            if d.used == True:
                ldihedrals.append(self.pd2ric(d.atoms))
        self.adihedrals = np.array(ldihedrals) + 1
        return

    def pd2ric(self, tor):
        ntor = 12 * [-1]
        l = [0,5,6,7]
        for i in range(4):
            ntor[l[i]] = tor[i]
        return ntor

    def get_weight_torsions(self, norm = True):
        self.wgh_dihedrals = []
        for d in self.pd.mol.dihedrals:
            weight = 0.0
            if (('torsion' in self.pd.mol.FF.variables.keys()) and
                    (d.type in self.pd.mol.FF.variables['torsion'].keys())):
                if norm == False:
                    weight = 1.0
                if norm == True:
                    weight = 1.0/len(self.pd.mol.FF.variables['torsion'][d.type][1])
            if d.used == True:
                self.wgh_dihedrals.append(weight)
        return
     
    def get_weights(self, norm = True):
        """
        Gets  the weights of the RICs

        :Parameters:
            - Norm : boolean
              If Norm = False all RICs inside the objective function get the same weight namely 1
              If Norm = True the weights per RIC type are normalized to one
        """
        ### bonds ###
        self.wgh_bonds = []
        for b in self.pd.mol.bonds:
            weight = 0.0
            if ((('bond' in self.pd.mol.FF.variables.keys()) and (b.type in self.pd.mol.FF.variables['bond'].keys())) or 
                    (('bondq' in self.pd.mol.FF.variables.keys()) and (b.type in self.pd.mol.FF.variables['bondq'].keys()))):
                if norm == False:
                    weight = 1.0
                if norm == True:
                    try:
                        weight = 1.0/len(self.pd.mol.FF.variables['bond'][b.type][1])
                    except KeyError:
                        weight = 1.0/len(self.pd.mol.FF.variables['bondq'][b.type][1])
            if b.used == True:
                self.wgh_bonds.append(weight)
        ### angles ###
        self.wgh_angles = []
        self.wgh_linangles = []
        for a in self.pd.mol.angles:
            weight = 0.0
            if ((('angle' in self.pd.mol.FF.variables.keys()) and (a.type in self.pd.mol.FF.variables['angle'].keys()))
                  or (('angleq' in self.pd.mol.FF.variables.keys()) and (a.type in self.pd.mol.FF.variables['angleq'].keys()))
                  or (('anglef' in self.pd.mol.FF.variables.keys()) and (a.type in self.pd.mol.FF.variables['anglef'].keys()))):
                if norm == False:
                    weight = 1.0
                if norm == True:
                    if (('anglef' in self.pd.mol.FF.variables.keys()) and (a.type in self.pd.mol.FF.variables['anglef'].keys())):
                        ###weight = 1.0/len(self.pd.mol.FF.variables['anglef'][a.type][1])
                        ###ADD
                        try:
                            weight = 1.0/len(self.pd.mol.FF.variables['anglef'][a.type][1])
                        except ZeroDivisionError:
                            raise ZeroDivisionError(
                            'Tip: either angle or anglef for the same triplet of atoms\n %s \n %s'
                            % (   str(self.pd.mol.FF.variables['anglef'][a.type][:]), str(a.type)   )
                            )
                        ###DDA
                    elif (('angleq' in self.pd.mol.FF.variables.keys()) and (a.type in self.pd.mol.FF.variables['angleq'].keys())):
                        weight = 1.0/len(self.pd.mol.FF.variables['angleq'][a.type][1])
                    else:
                        weight = 1.0/len(self.pd.mol.FF.variables['angle'][a.type][1])
            if a.used == True:
                if a.lin == False:  
                     self.wgh_angles.append(weight)
                if a.lin == True:
                     self.wgh_linangles.append(weight)
        ### oop ###
        self.wgh_oops = []
        for o in self.pd.mol.oops:
            weight = 0.0
            if (('opbend' in self.pd.mol.FF.variables.keys()) and (o.type in self.pd.mol.FF.variables['opbend'].keys())):
                if norm == False:
                    weight = 1.0
                if norm == True:
                    if o.tink == False:
                        #continue
                        weight = 1.0/(3*len(self.pd.mol.FF.variables['opbend'][o.type][1]))
                    else:
                        weight = 1.0/len(self.pd.mol.FF.variables['opbend'][o.type][1])
            self.wgh_oops.append(weight)
        ### strbends ###
        strbends = []
        self.wgh_strbends = []
        if 'strbnd' in self.pd.mol.FF.variables.keys():
            for i in range(len(self.pd.mol.FF.variables['strbnd'].keys())):
                 atomkey = self.pd.mol.FF.variables['strbnd'].keys()[i]
                 for a in self.pd.mol.angles:
                     if a.type == atomkey:
                         strbends.append(a.atoms)
                         if norm == True:
                             self.wgh_strbends.append(1.0/len(self.pd.mol.FF.variables['strbnd'][atomkey][1]))
                         else:
                             self.wgh_strbends.append(1.0)
        self.strbends = np.array(strbends) + 1 
        ### dihedrals  ###
        if self.fragtor == False:
            self.get_weight_torsions()
            return
        self.wgh_dihedrals = []
        for i in range(len(self.adihedrals)):
            weight = 0.0
            length = len(self.adihedrals[i])
            a1 = self.adihedrals[i][(length/2) -1]-1
            a2 = self.adihedrals[i][(length/2)]-1
            types = []
            types.append(self.pd.get_atomtypes()[a1])
            types.append(self.pd.get_atomtypes()[a2])
            if types in self.notwgh:
                continue
            count = 0
            for d in self.pd.mol.dihedrals:
                 if (((d.atoms[1] == a1) and (d.atoms[2] == a2))
                     or ((d.atoms[1] == a2) and (d.atoms[2] == a1))):
                     if (('torsion' in self.pd.mol.FF.variables.keys()) and (d.type in self.pd.mol.FF.variables['torsion'].keys())):
                         ### COUNT BONDS ###
                         if norm == False:
                             weight = 1.0
                         if norm == True:
                             for b in self.pd.mol.bonds:
                                 if (((self.pd.get_atomtypes()[b.atoms[0]] == types[0]) and (self.pd.get_atomtypes()[b.atoms[1]] == types[1])) or 
                                         ((self.pd.get_atomtypes()[b.atoms[0]] == types[1]) and (self.pd.get_atomtypes()[b.atoms[1]] == types[0]))):
                                     count +=1
                                    #weight = 1.0/len(self.pd.mol.FF.variables['torsion'][d.type][1])
                             weight = 1.0/count
                         break
            self.wgh_dihedrals.append(weight)
        return

    def check_lin(self, atoms, tresh = 3.0):
        """
        checks wether an angle between 3 atoms is linear

        :Parameters:
           - atoms:list
             list of indices of the involved atoms
        """
        #phi = tools.geometry(self.pd.get_xyz()).get_angle(atoms)
        phi = tools.geometry(self.refxyz).get_angle(atoms)
        if abs(180.0 - phi) < tresh:
            return True
        else:
            return False

    def choose_lin(self, atoms):
        """
        atoms (tuple) : triplet of atoms belonging to the linear angle
        known problems: it can assign rotating atoms wtr. the triplet (e.g. H_Me)
        and could cause instability (not yet tested)
        xyz: atoms coordinates
        vec: vector along the external angles
        xyzmod
        """
        xyz = self.pd.mol.xyz
        natoms = self.pd.mol.natoms
        norm = np.linalg.norm
        vec = [xyz[atoms[2]][i] - xyz[atoms[0]][i] for i in xrange(3)]
        xyzmod = [vec if ixyz in atoms else exyz for ixyz,exyz in enumerate(xyz)]
        #crsxyz = np.cross(vec, xyz)
        crsxyz = np.cross(vec, xyzmod)
        sinang2 = [None] * natoms
        for i in xrange(natoms):
            if norm(crsxyz[i]) == 0.:
                sinang2[i] = np.inf
            else:
                sinang2[i] = norm(xyz[i])*norm(vec)/norm(crsxyz[i])
        atomref = np.argmin(sinang2) + 1 ###due to inconsitency
        return atomref

    def write_ric(self):
        """
        Gets ric form pydlpoly and writes them to ric.py
        """
        for i in range(np.shape(self.abonds)[0]):
            self.pd.ric.add_stretch(self.abonds[i,:])
        for i in range(np.shape(self.aangles)[0]):
            self.pd.ric.add_in_bend(self.aangles[i,:])
        for i in range(np.shape(self.alinangles)[0]):
            self.pd.ric.add_lin_bend(self.alinangles[i,:],
                    self.lin_ref_atoms[i])
        for i in range(np.shape(self.aoops)[0]):
            self.pd.ric.add_out_bend(self.aoops[i,:])
        for i in range(len(self.adihedrals)):
            self.pd.ric.add_torsion(self.adihedrals[i][:])
        self.pd.ric.add_eckart()
        print '##### defined rics #####'
        for i in range(len(self.pd.ric._stretches)):
            print 'str', self.pd.ric._stretches[i]
        for i in range(len(self.pd.ric._in_bends)):
            print 'ibe', self.pd.ric._in_bends[i]
        for i in range(len(self.pd.ric._out_bends)):
            print 'obe', self.pd.ric._out_bends[i]
        for i in range(len(self.pd.ric._lin_bends)):
            print 'lin', self.pd.ric._lin_bends[i]
        for i in range(len(self.pd.ric._torsions)):
            print 'tor', self.pd.ric._torsions[i]
        print '########################'
        return

    def get_formatted_torsions(self):
        unformatted = self.pd.ric._torsions
        formatted = []
        for i in unformatted:
            formatted.append(filter(lambda a: a !=0, i))
        return formatted

