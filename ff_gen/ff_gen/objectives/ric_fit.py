####################################################################
#
#  OUTDATED! This module has been moved to ric_fit2 and is not
#  longer updated. Please use ric_fit2
#
###################################################################


"""
This file implements a ric_fit class.
It is inheriting all features from RedIntCoords in ric.py
and adds features to load a respective reference structure and Hessian.
In addition, a weight matrix is held to evaluate various kinds of weighted mean
square deviations to be used as ingredients to fittness values
"""

import string
import numpy as np
from ric import RedIntCoords 
import copy
import os
import turbomole
import refclass
import time
import sys

ricdic = {'str':[2,0],
        'ibe' : [3,1],
        'obe' : [4,2],
        'lbe' : [3,3],
        'tor' : [12,4]}


def revlist(l):
    nl = copy.copy(l)
    nl.revert()
    return nl

class ric_fit:
    """
    class to compute redundant internal coordinates (ric)
    by using the inherited ric module and to compute deviations from a corresponding
    reference.
    """
    
    def __init__(self, tag = 'primary', disp = 0.001, verbose = False, fragtor=False, lin_ref = None, notwgh=[],
            oweight = [1,1,1], absolute = True):
        self.tag = tag
        self.disp = disp
        self.verbose = verbose
        self.notwgh = notwgh
        if lin_ref != None:
            self.read_lin_ref(lin_ref)
        self.oweight = np.array(oweight)
        self.fragtor = fragtor
        self.absolute = absolute
        sys.stderr.write('OUTDATED! Please use ric_fit2. This module is not longer updated\n')
        return
    
    def initialize(self, pd, ref):
        self.pd = pd
        self.ref =ref
        self.get_ric()
        # init ric itself
        if hasattr(self.pd, 'ric') == False:
            self.pd.ric = RedIntCoords() 
            self.write_ric()
            self.pd.ric.setup(self.pd.get_masses())
        self.setup()
        self.get_weights(norm = True)
        self.write_weight()
#        self.ixyz = copy.deepcopy(self.pd.get_xyz())
        self.ixyz = copy.deepcopy(self.ref(info = 'coord', branch = 'hessians', tag = self.tag))
        return

    def setup(self):
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
        # Allocate a reference structure/hessian array
        self.ref_val_str = np.zeros([self.nstretch],dtype="float64")
        self.wgt_val_str = np.zeros([self.nstretch],dtype="float64")
        self.ref_val_ibe = np.zeros([self.nin_bend],dtype="float64")
        self.wgt_val_ibe = np.zeros([self.nin_bend],dtype="float64")
        self.ref_val_obe = np.zeros([self.nout_bend],dtype="float64")
        self.wgt_val_obe = np.zeros([self.nout_bend],dtype="float64")
        self.ref_val_lbe = np.zeros([self.nlin_bend],dtype="float64")
        self.wgt_val_lbe = np.zeros([self.nlin_bend],dtype="float64")
        self.ref_val_tor = np.zeros([self.ntorsion],dtype="float64")
        self.wgt_val_tor = np.zeros([self.ntorsion],dtype="float64")
        self.ref_hes = np.zeros([self.pd.ric.nric, self.pd.ric.nric], dtype="float64")
        self.wgt_hes = np.zeros([self.pd.ric.nric, self.pd.ric.nric], dtype="float64")
        # set some defaults
        self.fact_str = 1.0
        self.fact_ibe = 1.0
        self.fact_obe = 1.0
        self.fact_lbe = 1.0
        self.fact_tor = 1.0
        self.cycle = 0
        self.fdiagnostics = open('ric_fit_%s.punch' % self.tag, 'w')
        # We use here the original convention to compare Hessian elements in milidyn*A
        self.fact_hes = 1.0/143.88
        # generate referece data        
        self.generate_reference(self.pd.get_cell())
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
#        assert self._setup
        if ric_type=="str":
            assert len(ind) == 2
            try:
                i = self.pd.ric._stretches.index(ind)
            except ValueError:
                ind.reverse()
                try:
                    i = self.pd.ric._stretches.index(ind)
                except ValueError:
                    return -1, -1
            iglob = self.first_str+i
            return i, iglob
        elif ric_type=="ibe":
            assert len(ind) == 3
            try:
                i = self.pd.ric._in_bends.index(ind)
            except ValueError:
                ind.reverse()
                try:
                    i = self.pd.ric._in_bends.index(ind)
                except ValueError:
                    return -1, -1
            iglob = self.first_ibe+i
            return i, iglob
        elif ric_type=="obe":
            assert len(ind) == 4
            try:
                i = self.pd.ric._out_bends.index(ind)
            except ValueError:
                ind = ind[:2]+revlist(ind[2:])
                try:
                    i = self.pd.ric._out_bends.index(ind)
                except ValueError:
                    return -1, -1
            iglob = self.first_obe+i
            return i, iglob
        elif ric_type=="lbe":
            assert len(ind) == 3
            try:
                i = self.pd.ric._lin_bends.index(ind)
            except ValueError:
                ind.reverse()
                try:
                    i = self.pd.ric._lin_bends.index(ind)
                except ValueError:
                    return -1, -1
            iglob = self.first_lbe+i
            return i, iglob
        elif ric_type=="tor":
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
        else:
            raise "ValueError", "Unknown ric type %s" % ric_type
        return

    def map_ric2(self, ric_type, ind):
        assert len(ind) == ricdic[ric_type][0]
        try:
            i = self.riclist[ricdic[ric_type][1]].index(ind)
#            i = self.pd.ric._lin_bends.index(ind)
        except ValueError:
            ind.reverse()
            try:
                i = self.riclist[ricdic[ric_type][1]].index(ind)
            except ValueError:
                return -1, -1
        iglob = self.firstrics[ricdic[ric_type][1]]+i
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
  
    def get_ric(self):
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
        for a in self.pd.mol.angles:
            if a.used == True:
                if self.check_lin(a.atoms) == False:
                    langles.append(a.atoms)
                    a.lin = False
                if self.check_lin(a.atoms) == True:
                    if a.atoms in self.lin_dict.keys():
                        self.lin_ref_atoms.append(self.lin_dict[a.atoms])
                        llinangles.append(a.atoms)
                        a.lin = True
                    else:
                        print 'no ref atom definition found for lin bend %d' % a.atoms[1]+1
                        raise IOError
        self.aangles = np.array(langles) + 1
        self.alinangles = np.array(llinangles) + 1
        ### oop ###
        loops = []
        i = 0
        for o in self.pd.mol.oops:
#      if o.used == True:
#            if (i/3.0 - i/3) == 0.0:
#                a = o.atoms[0]
#                b = o.atoms[1]
#                c = o.atoms[2]
#                d = o.atoms[3]
#                if self.check_lin([c,a,d]) == False:
#                    loops.append([a,b,c,d])
#                if self.check_lin([b,a,d]) == False:
#                    loops.append([a,c,b,d])
#                if self.check_lin([b,a,c]) == False:
#                    loops.append([a,d,b,c])
#            i += 1
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
                        weight = 1.0/len(self.pd.mol.FF.variables['anglef'][a.type][1])
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
        ### strbends ###
#        strbends = []
#        self.wgh_strbends = []
#        if 'strbnd' in self.pd.mol.FF.variables.keys():
#            for i in range(len(self.pd.mol.FF.variables['strbnd'].keys())):
#                 atomkey = self.pd.mol.FF.variables['strbnd'].keys()[i]
#                 for a in self.pd.mol.angles:
#                     if a.type == atomkey:
#                         strbends.append(a.atoms)
#                         if norm == True:
#                             self.wgh_strbends.append(1.0/len(self.pd.mol.FF.variables['strbnd'][atomkey][1]))
#                         else:
#                             self.wgh_strbends.append(1.0)
#        self.strbends = np.array(strbends) + 1 
        return

    def check_lin(self, atoms, tresh = 0.0038):
        """
        checks wether an angle between 3 atoms is linear

        :Parameters:
           - atoms:list
             list of indices of the involved atoms
        """
        apex_1 = self.pd.get_xyz()[atoms[0],:]
        apex_2 = self.pd.get_xyz()[atoms[2],:]
        central = self.pd.get_xyz()[atoms[1],:]
        r1 = apex_1 - central
        r2 = apex_2 - central
        nr1 = r1/np.sqrt((r1*r1).sum())
        nr2 = r2/np.sqrt((r2*r2).sum())
        s = np.dot(nr1, nr2)
        if 1+s < tresh:
        #if 1+s < 0.01:
            return True
        else:
            return False

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

    def write_weight(self):
        """
        Gets the weights corresponing to the rics and writes them to ric.py
        """
        for j in range(len(self.wgh_bonds)):
            weight = self.wgh_bonds[j]
            i, iglob = self.map_ric("str", list(self.abonds[j,:]))
            self.wgt_val_str[i]        = weight
            self.wgt_hes[iglob, iglob] = weight
        for j in range(len(self.wgh_angles)):
            weight = self.wgh_angles[j]
            i, iglob = self.map_ric("ibe", list(self.aangles[j,:]))
            self.wgt_val_ibe[i]        = weight
            self.wgt_hes[iglob, iglob] = weight
        for j in range(len(self.wgh_linangles)):
            weight = self.wgh_linangles[j]
            i, iglob = self.map_ric("lbe", list(self.alinangles[j,:]))
#            self.wgt_val_ibe[i]        = weight
            self.wgt_hes[iglob, iglob] = weight/2.0
            self.wgt_hes[iglob+1, iglob+1] = weight/2.0
        for j in range(len(self.wgh_oops)):
            weight = self.wgh_oops[j]
            i, iglob = self.map_ric("obe", list(self.aoops[j,:]))
            self.wgt_val_obe[i]        = weight
            self.wgt_hes[iglob, iglob] = weight
        for j in range(len(self.wgh_dihedrals)):
            if self.fragtor == False:
                i, iglob = self.map_ric2("tor", list(self.adihedrals[j,:]))
            else:
                ind = self.adihedrals[j]
                na = len(self.adihedrals[j])
                ind = ind[na/2-1:na/2+1]
                i, iglob = self.map_ric("tor", list(ind))
            weight = self.wgh_dihedrals[j]
            self.wgt_hes[iglob, iglob] = weight
        #self.wgt_hes[:-6,:-6] = 1
        ### str - bnd ###
        for j in range(np.shape(self.strbends)[0]):
            weight = self.wgh_strbends[j]/6.0
            i, iglob_bend = self.map_ric("ibe", list(self.strbends[j,:]))
            i, iglob_str1 = self.map_ric("str", list(self.strbends[j,:2]))
            i, iglob_str2 = self.map_ric("str", list(self.strbends[j,1:]))
            self.wgt_hes[iglob_bend, iglob_str1] += weight
            self.wgt_hes[iglob_str1, iglob_bend] += weight
            self.wgt_hes[iglob_bend, iglob_str2] += weight
            self.wgt_hes[iglob_str2, iglob_bend] += weight
            self.wgt_hes[iglob_str1, iglob_str2] += weight
            self.wgt_hes[iglob_str2, iglob_str1] += weight
        return

    def generate_reference(self, cell):
#        self.pd.ric.construct_b_matrix(cell, self.ref('coord'))
        self.pd.ric.construct_b_matrix(cell, self.ref(info = 'coord', branch = 'hessians', tag = self.tag))
        invb, rankb = self.pd.ric.invert_b_matrix()
        if rankb < self.ref(info = 'natoms', branch = 'system'):
            print 'ERROR: Rank is too small'
            raise ValueError
        self.pd.ric.project_hessian(self.ref(info = 'hessian', branch = 'hessians', tag = self.tag))
        self.ref_val_str = copy.deepcopy(self.pd.ric.get_val_stretches())
        self.ref_val_ibe = copy.deepcopy(self.pd.ric.get_val_in_bends())
        self.ref_hes     = copy.deepcopy(self.pd.ric.get_ric_hessian())
        return

    def read_ric(self, filename):
        """
        Reads ric from the original "new_reference.dat" file from ff_generator
    
        the rics must be ordered in the canonical form (stretch, bend, wag, torsion)
        This function must be called before setup.
    
        :Parameters:
    
        - filename:    Name of the input file (usually "new_reference.dat")
    
        """
        f = open(filename, "r")
        line = f.readline()
        if line[:3] != "###":
            raise IOError, "This does not look like a traditional reference file"
        line = f.readline()
        stop = False
        while not stop:
            sline = line.split()
            if sline[1] == "stretch":
                ind = map(string.atoi, sline[2:4])
                self.add_stretch(ind)
            elif sline[1] == "bend":
                ind = map(string.atoi, sline[2:5])
                self.add_in_bend(ind)
            elif sline[1] == "wag":
                ind = map(string.atoi, sline[2:6])
                self.add_out_bend(ind)
            elif sline[1] == "torsion":
                for i in xrange(2, len(sline)):
                    if sline[i].count(".") == 1: break
                ind = map(string.atoi, sline[2:i])
                self.add_torsion(ind)
            else:
                raise IOError, "Unknown RIC type in reference file"
            line = f.readline()
            if len(line)==0: stop=True
            # in the refernece_dat file the Hessian elements start with a line with "###" -> stop here
            if line[:3] == "###": stop=True
        f.close()
        # HACK
        # assumption ... non peridic .. fix
        self.add_eckart_rots()
        return
                            
                            
    def read_ref_from_refdat(self, filename, conv_angle=np.pi/180.0, conv_energy = 143.88, conv_length=1.0, read_fit_flag=True):
        """
        Reads reference data from the original "new_reference.dat" file from ff_generator
    
        Currently torsion and out_bends do not support reading a value in
        Off-diagonal values of the Hessian are not supported
        Must be called after setup
    
        :Parameters:
    
        - filename:         Name of the input file (usually "new_reference.dat")
        - [conv_angle]:     Conversion for angles to give radians, defualt is to read degree values
        - [conv_energy]:    Conversion to give the (pydlpoly default) kcal/mol from the old default milidyne*A
        - [conv_length]:    Conversion factor to Angstrom .. default is 1.0
        - [read_fit_flag]:  If True (default) read in fit flag and set weight to 1.0 for this RIC
    
        """
        assert self._setup
        #
        f = open(filename, "r")
        line = f.readline()
        if line[:3] != "###":
            raise IOError, "This does not look like a traditional reference file"
        line = f.readline()
        stop = False
        while not stop:
            sline = line.split()
            if sline[1] == "stretch":
                ind = map(string.atoi, sline[2:4])
                ref_val = string.atof(sline[4])*conv_length
                ref_2nd = string.atof(sline[5])*conv_energy
                weight = 0.0
                if len(sline) > 6:
                    if sline[6] == "fit": weight = 1.0
                i, iglob = self.map_ric("str", ind)
                if i >=0:
                    self.ref_val_str[i]        = ref_val
                    self.ref_hes[iglob, iglob] = ref_2nd
                    print "read RIC stretch %20s ref values %10.5f %10.5f" % (ind, ref_val, ref_2nd) 
                    if read_fit_flag:
                        self.wgt_val_str[i]        = weight
                        self.wgt_hes[iglob, iglob] = weight
                        if weight >0: print "        RIC is fitted"
                else:
                    print "RIC stretch %s is not defiend" % ind
            elif sline[1] == "bend":
                ind = map(string.atoi, sline[2:5])
                ref_val = string.atof(sline[5])*conv_angle
                ref_2nd = string.atof(sline[6])*conv_energy
                weight = 0.0
                if len(sline) > 7:
                    if sline[7] == "fit": weight = 1.0
                i, iglob = self.map_ric("ibe", ind)
                if i >=0:
                    self.ref_val_ibe[i]        = ref_val
                    self.ref_hes[iglob, iglob] = ref_2nd
                    print "read RIC bend    %20s ref values %10.5f %10.5f" % (ind, ref_val, ref_2nd) 
                    if read_fit_flag:
                        self.wgt_val_ibe[i]        = weight
                        self.wgt_hes[iglob, iglob] = weight
                        if weight >0: print "        RIC is fitted"
                else:
                    print "RIC bend %s is not defiend" % ind
            elif sline[1] == "wag":
                ind = map(string.atoi, sline[2:6])
                ref_2nd = string.atof(sline[6])*conv_energy
                weight = 0.0
                if len(sline) > 7:
                    if sline[7] == "fit": weight = 1.0
                i, iglob = self.map_ric("obe", ind)
                if i >=0:
                    self.ref_hes[iglob, iglob] = ref_2nd
                    print "read RIC wag     %20s ref values %10.5f" % (ind, ref_2nd) 
                    if read_fit_flag:
                        self.wgt_hes[iglob, iglob] = weight
                        if weight >0: print "        RIC is fitted"
                else:
                    print "RIC wag %s is not defiend" % ind
            elif sline[1] == "torsion":
                for i in xrange(2, len(sline)):
                    if sline[i].count(".") == 1: break
                ind = map(string.atoi, sline[2:i])
                na = len(ind)
                ind = ind[na/2-1:na/2+1]
                weight = 0.0
                if len(sline) > i+1:
                    if sline[i+1] == "fit": weight = 1.0
                ref_2nd = string.atof(sline[i])*conv_energy
                i, iglob = self.map_ric("tor", ind)
                if i >=0:
                    self.ref_hes[iglob, iglob] = ref_2nd
                    print "read RIC torsion %20s ref values %10.5f" % (ind, ref_2nd) 
                    if read_fit_flag:
                        self.wgt_hes[iglob, iglob] = weight
                        if weight >0: print "        RIC is fitted"
                else:
                    print "RIC torsion %s is not defiend" % ind
            else:
                raise IOError, "Unknown RIC type in reference file"
            line = f.readline()
            if len(line)==0: stop=True
            # in the refernece_dat file the Hessian elements start with a line with "###" -> stop here
            if line[:3] == "###": stop=True
        f.close()
        return

    def __call__(self):
        call_start = time.clock()
        self.pd.set_xyz(self.ixyz)
        self.pd.set_atoms_moved()
        self.pd.MIN_lbfgs(0.001)
        #self.pd.write_tinker_xyz('actual.xyz')
        const_start = time.clock()
        self.pd.ric.construct_b_matrix(self.pd.get_cell(), self.pd.get_xyz())
        const_ende = time.clock()
        invert_start = time.clock()
        self.pd.ric.invert_b_matrix()
        invert_ende = time.clock()
        hess_start = time.clock()
        hessian = self.pd.calc_hessian(disp = self.disp)
        hess_ende = time.clock()
        #self.pd.ric.project_hessian(self.pd.calc_hessian(disp = self.disp))
        project_start = time.clock()
        self.pd.ric.project_hessian(hessian)
        project_ende = time.clock()
        self.calc_msd()
#        self.msd = self.msd_str + self.msd_ibe + self.msd_hes
        self.msd = np.sum(self.oweight*np.array([self.msd_str,self.msd_ibe,self.msd_hes]))
        self.fdiagnostics.write("%s %6.6f %6.6f %6.6f % 6.6f\n" % (self.cycle, self.msd, 
            self.msd_str, self.msd_ibe, self.msd_hes))
        if self.verbose == True:
            self.print_diagonals()
        self.cycle += 1
        call_ende = time.clock()
#        print 'time gesamt: ', call_ende-call_start
#        print 'time constr: ', const_ende-const_start
#        print 'time invert: ', invert_ende-invert_start
#        print 'time hess:   ', hess_ende-hess_start
#        print 'time projec: ', project_ende-project_start
        return self.msd, [[self.msd, self.msd_str, self.msd_ibe, self.msd_hes]]


    def calc_msd(self):
        """ 
        Calculates the mean square deviations and returns them as a numpy array
        """
        # stretch
        if self.wgt_val_str.sum() > 0:
            delt = (self.pd.ric.get_val_stretches()-self.ref_val_str)*self.fact_str
            wdelt = delt*self.wgt_val_str
            self.msd_str = (wdelt*wdelt).sum()/self.wgt_val_str.sum()
            self.max_str = wdelt.max()
            self.amsd_str = delt * delt
        else:
            self.msd_str = 0.0
            self.amsd_str = 0.0
        # bend
        if self.wgt_val_ibe.sum() > 0:
            delt = (self.pd.ric.get_val_in_bends()-self.ref_val_ibe)*self.fact_ibe
            wdelt = delt*self.wgt_val_ibe
            self.msd_ibe = (wdelt*wdelt).sum()/self.wgt_val_ibe.sum()
            self.max_ibe = wdelt.max()
            self.amsd_ibe = delt * delt
        else:
            self.msd_ibe = 0.0
            self.amsd_ibe = 0.0
        # wag ... not compared
        # torsion ... to be implemented
        # hessian (currently all hessian elements are compared
        if self.absolute:
            delt = (self.pd.ric.get_ric_hessian()-self.ref_hes)*self.fact_hes
        else:
            delt = abs(1-(self.pd.ric.get_ric_hessian()/self.ref_hes))
        wdelt = delt*self.wgt_hes
        self.msd_hes = (wdelt*wdelt).sum()/self.wgt_hes.sum()
        self.max_hes = wdelt.max()
        self.amsd_hes = delt * delt
        return
    
    def print_diagonals(self):
        rics = self.pd.ric._stretches + self.pd.ric._in_bends + \
                self.pd.ric._out_bends +self.pd.ric._lin_bends + self.pd.ric._torsions +\
                ['e'] +['e'] +['e'] +['e'] +['e'] +['e']
        for i in range(len(rics)):
            diff = self.ref_hes[i, i]-self.pd.ric.get_ric_hessian()[i,i]
            ratio = abs(1-(self.pd.ric.get_ric_hessian()[i,i]/self.ref_hes[i,i]))*100
            print i, rics[i], self.ref_hes[i, i], self.pd.ric.get_ric_hessian()[i,i], \
                    diff, ratio, self.wgt_hes[i,i]
        return
 
