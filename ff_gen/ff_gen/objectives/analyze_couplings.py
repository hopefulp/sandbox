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

def revlist(l):
    nl = copy.copy(l)
    nl.revert()
    return nl

class analyze_couplings:
    """
    class to compute redundant internal coordinates (ric)
    by using the inherited ric module and to compute deviations from a corresponding
    reference.
    """
    
    def __init__(self, tag = 'primary'):
        self.tag = tag
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
        # We use here the original convention to compare Hessian elements in milidyn*A
        self.fact_hes = 1.0/143.88
        # generate referece data        
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
        for a in self.pd.mol.angles:
            if a.used == True:
                if self.check_lin(a.atoms) == False:
                    langles.append(a.atoms)
                    a.lin = False
                if self.check_lin(a.atoms) == True:
                    llinangles.append(a.atoms)
                    a.lin = True
        self.aangles = np.array(langles) + 1
        self.alinangles = np.array(llinangles) + 1
        ### oop ###
        loops = []
        i = 0
        for o in self.pd.mol.oops:
#      if o.used == True:
            if (i/3.0 - i/3) == 0.0:
                a = o.atoms[0]
                b = o.atoms[1]
                c = o.atoms[2]
                d = o.atoms[3]
                loops.append([a,b,c,d])
                loops.append([a,c,b,d])
                loops.append([a,d,b,c])
            i += 1
        self.aoops = np.array(loops) +1
        ### dihedrals ###
        ldihedrals = []
        for i in range(len(lbonds)):
            a1 = lbonds[i][0]
            a2 = lbonds[i][1]
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
     

    def check_lin(self, atoms):
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
        if 1+s < 0.0038:
            return True
        else:
            return False

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
                  or (('anglef' in self.pd.mol.FF.variables.keys()) and (a.type in self.pd.mol.FF.variables['anglef'].keys()))):
                if norm == False:
                    weight = 1.0
                if norm == True:
                    if (('anglef' in self.pd.mol.FF.variables.keys()) and (a.type in self.pd.mol.FF.variables['anglef'].keys())):
                        weight = 1.0/len(self.pd.mol.FF.variables['anglef'][a.type][1])
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
                    weight = 1.0/(3*len(self.pd.mol.FF.variables['opbend'][o.type][1]))
            self.wgh_oops.append(weight)
        ### dihedrals  ###
        self.wgh_dihedrals = []
        for i in range(len(self.adihedrals)):
            weight = 0.0
            length = len(self.adihedrals[i])
            a1 = self.adihedrals[i][(length/2) -1]-1
            a2 = self.adihedrals[i][(length/2)]-1
            types = []
            types.append(self.pd.get_atomtypes()[a1])
            types.append(self.pd.get_atomtypes()[a2])
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
        return



    def write_ric(self):
        """
        Gets ric form pydlpoly and writes them to ric.py
        """
        for i in range(np.shape(self.abonds)[0]):
            self.pd.ric.add_stretch(self.abonds[i,:])
        for i in range(np.shape(self.aangles)[0]):
            self.pd.ric.add_in_bend(self.aangles[i,:])
        for i in range(np.shape(self.alinangles)[0]):
            self.pd.ric.add_lin_bend(self.anlinangles[i,:])
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
            self.wgt_hes[iglob, iglob] = weight
        for j in range(len(self.wgh_oops)):
            weight = self.wgh_oops[j]
            i, iglob = self.map_ric("obe", list(self.aoops[j,:]))
#            self.wgt_val_obe[i]        = weight
            self.wgt_hes[iglob, iglob] = weight
        for j in range(len(self.wgh_dihedrals)):
            ind = self.adihedrals[j]
            na = len(self.adihedrals[j])
            ind = ind[na/2-1:na/2+1]
            weight = self.wgh_dihedrals[j]
            i, iglob = self.map_ric("tor", ind)
#            self.wgt_val_obe[i]        = weight
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


    def __call__(self, rictype, ric):
#        self.pd.ric.construct_b_matrix(cell, self.ref('coord'))
        self.pd.ric.construct_b_matrix(self.pd.get_cell(), self.ref(info = 'coord', branch = 'hessians', tag = self.tag))
        invb, rankb = self.pd.ric.invert_b_matrix()
        if rankb < self.ref(info = 'natoms', branch = 'system'):
            print 'ERROR: Rank is too small'
            raise ValueError
        self.pd.ric.project_hessian(self.ref(info = 'hessian', branch = 'hessians', tag = self.tag))
        self.ric_hes     = copy.deepcopy(self.pd.ric.get_ric_hessian())
        self.pd.ric.construct_b_matrix(self.pd.get_cell(), self.pd.get_xyz())
        invb, rankb = self.pd.ric.invert_b_matrix()
        self.pd.ric.project_hessian(self.pd.calc_hessian())
        i, iglob = self.map_ric(rictype, ric)
        rics = self.pd.ric._stretches + self.pd.ric._in_bends + \
                self.pd.ric._out_bends + self.pd.ric._torsions +\
                ['e'] +['e'] +['e'] +['e'] +['e'] +['e']
        print rictype, ric
        for i in range(len(rics)):
            diff =self.ric_hes[i, iglob]-self.pd.ric.get_ric_hessian()[i,iglob]
            #if abs(diff) > 0.5:
            if abs(diff) > 1.0:
                print i, rics[i], self.ric_hes[i, iglob], self.pd.ric.get_ric_hessian()[i,iglob], \
                        diff, self.wgt_hes[i,iglob]
        return


    def calc_msd(self):
        """ 
        Calculates the mean square deviations and returns them as a numpy array
        """
#        assert self._setup
#        assert self._eval
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
        delt = (self.pd.ric.get_ric_hessian()-self.ref_hes)*self.fact_hes
        wdelt = delt*self.wgt_hes
        self.msd_hes = (wdelt*wdelt).sum()/self.wgt_hes.sum()
        self.max_hes = wdelt.max()
        self.amsd_hes = delt * delt
        return
            
    
