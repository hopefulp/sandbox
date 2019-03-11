"""
This file implements a ric_fit class.
It is inheriting all features from RedIntCoords in ric.py
and adds features to load a respective reference structure and Hessian.
In addition, a weight matrix is held to evaluate various kinds of weighted mean
square deviations to be used as ingredients to fittness values
"""

import string
import numpy as np
from scipy import signal
import copy
import os
import time 

from ff_gen.ric import RedIntCoords 
import ff_gen.turbomole as turbomole
import ff_gen.refclass as refclass
import  pydlric
import hessutils

rad2deg = 180.0/np.pi

def revlist(l):
    nl = copy.copy(l)
    nl.revert()
    return nl

class ric_fit(pydlric.pydlric):
    """
    class to compute redundant internal coordinates (ric)
    by using the inherited ric module and to compute deviations from a corresponding
    reference.
    """
    
    def __init__(self, tag = 'primary', disp = 0.001, verbose = False, fragtor=False, lin_ref = None, notwgh=[],
            oweight = [1,1,0,1], absolute = True):
        self.tag = tag
        self.disp = disp
        self.verbose = verbose
        self.notwgh = notwgh
        if lin_ref != None:
            self.read_lin_ref(lin_ref)
        self.oweight = np.array(oweight)
        self.fragtor = fragtor
        self.absolute = absolute
        return
    
    def initialize(self, pd, ref):
        self.pd = pd
        self.ref =ref
        self.refxyz = self.ref(info = 'coord', branch = 'hessians', tag = self.tag)
        self.get_ric()
        # init ric itself
        if hasattr(self.pd, 'ric') == False:
            self.pd.ric = RedIntCoords() 
            self.write_ric()
            self.pd.ric.setup(self.pd.get_masses())
        self.prepare()
        self.setup()
        self.get_weights(norm = True)
        self.write_weight()
        self.ixyz = copy.deepcopy(self.ref(info = 'coord', branch = 'hessians', tag = self.tag))
        self.hess=hessutils.hessian(self.pd)
        return

    def setup(self):
        """
        Setup the RICs 
        
        Call after adding all rics
        
        :Parameters:
        
            - masses : numpy array [natoms] defining the atomic masses
            
        """
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
        self.fact_tor = 180.0/np.pi
        self.cycle = 0
        self.fdiagnostics = open('ric_fit_%s.punch' % self.tag, 'w')
        # We use here the original convention to compare Hessian elements in milidyn*A
        self.fact_hes = 1.0/143.88
        # generate referece data        
        self.generate_reference(self.pd.get_cell())
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
            self.wgt_hes[iglob, iglob] = weight/2.0
            self.wgt_hes[iglob+1, iglob+1] = weight/2.0
        for j in range(len(self.wgh_oops)):
            weight = self.wgh_oops[j]
            i, iglob = self.map_ric("obe", list(self.aoops[j,:]))
            self.wgt_val_obe[i]        = weight
            self.wgt_hes[iglob, iglob] = weight
        for j in range(len(self.wgh_dihedrals)):
            if self.fragtor == False:
                i, iglob = self.map_ric("tor", list(self.adihedrals[j,:]))
            else:
                ind = self.adihedrals[j]
                na = len(self.adihedrals[j])
                ind = ind[na/2-1:na/2+1]
                i, iglob = self.map_fractor_ric(list(ind))
            weight = self.wgh_dihedrals[j]
            self.wgt_val_tor[i] = weight
            self.wgt_hes[iglob, iglob] = weight
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
        self.pd.ric.construct_b_matrix(cell, self.ref(info = 'coord', branch = 'hessians', tag = self.tag))
        invb, rankb = self.pd.ric.invert_b_matrix()
        if rankb < self.ref(info = 'natoms', branch = 'system'):
            print 'ERROR: Rank is too small'
            raise ValueError
        self.pd.ric.project_hessian(self.ref(info = 'hessian', branch = 'hessians', tag = self.tag))
        self.ref_val_str = copy.deepcopy(self.pd.ric.get_val_stretches())
        self.ref_val_ibe = copy.deepcopy(self.pd.ric.get_val_in_bends())
        self.ref_val_tor = copy.deepcopy(self.pd.ric.get_val_torsions())
        self.ref_hes     = copy.deepcopy(self.pd.ric.get_ric_hessian())
        return

    def __call__(self):
        call_start = time.clock()
        self.pd.set_xyz(self.ixyz)
        self.pd.set_atoms_moved()
        self.pd.MIN_lbfgs(0.001, verbose = False)
        #self.pd.write_tinker_xyz('actual.xyz')
        const_start = time.clock()
        self.pd.ric.construct_b_matrix(self.pd.get_cell(), self.pd.get_xyz())
        const_ende = time.clock()
        invert_start = time.clock()
        self.pd.ric.invert_b_matrix()
        invert_ende = time.clock()
        hess_start = time.clock()
#        hessian = self.pd.calc_hessian(disp = self.disp)
        hessian = self.hess(disp = self.disp)
        hess_ende = time.clock()
        project_start = time.clock()
        self.pd.ric.project_hessian(hessian)
        project_ende = time.clock()
        self.calc_msd()
#        self.msd = self.msd_str + self.msd_ibe + self.msd_hes
        self.msd = np.sum(self.oweight*np.array([self.msd_str,self.msd_ibe,self.msd_tor,self.msd_hes]))
#        self.fdiagnostics.write("%s %6.6f %6.6f %6.6f % 6.6f\n" % (self.cycle, self.msd, 
#            self.msd_str, self.msd_ibe, self.msd_hes))
        if self.verbose == True:
            self.print_diagonals()
            self.print_geometry()
        self.cycle += 1
        call_ende = time.clock()
#        print 'time gesamt: ', call_ende-call_start
#        print 'time constr: ', const_ende-const_start
#        print 'time invert: ', invert_ende-invert_start
#        print 'time hess:   ', hess_ende-hess_start
#        print 'time projec: ', project_ende-project_start
        return self.msd, [[self.msd, self.msd_str, self.msd_ibe, self.msd_tor, self.msd_hes]]


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
        # torsion ... to be fully implemented ### DO NOT REMOVE COMMENTS
        if self.wgt_val_tor.sum() > 0:
            delt = (self.pd.ric.get_val_torsions()-self.ref_val_tor)*self.fact_tor
            delt = 90*(signal.sawtooth(delt/rad2deg,0.5)+1) ###only multiplicity 1
            wdelt = delt*self.wgt_val_tor
            self.msd_tor = (wdelt*wdelt).sum()/self.wgt_val_tor.sum()
            self.max_tor = wdelt.max()
            self.amsd_tor = delt * delt
        else:
            self.msd_tor = 0.0
            self.amsd_tor = 0.0
        # wag ... not compared
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
                self.pd.ric._out_bends +self.pd.ric._lin_bends + self.get_formatted_torsions() +\
                ['e'] +['e'] +['e'] +['e'] +['e'] +['e']
        types = self.nstretch*['str']+self.nin_bend*['ibe']+self.nout_bend*['obe']+self.nlin_bend*['lin']\
                + self.ntorsion*['tor'] + 6*['eck']
        print "# %6s %5s %20s %12s %12s %12s %12s %12s" %('idx', 'type', 'atoms', 'ref', 'ff', 'abs diff', 'rel diff', 'weight')
        for i in range(len(rics)):
            diff = self.ref_hes[i, i]-self.pd.ric.get_ric_hessian()[i,i]
            ratio = abs(1-(self.pd.ric.get_ric_hessian()[i,i]/self.ref_hes[i,i]))*100
            print "%6d %5s %20s %12.4f %12.4f %12.4f %12.4f %12.4f" %\
                    (i, types[i], str(rics[i]), self.ref_hes[i, i], self.pd.ric.get_ric_hessian()[i,i], \
                    diff, ratio, self.wgt_hes[i,i])
        return
    
    def print_geometry(self):
        print '---------stretches----------'
        print "# %6s %20s %12s %12s %12s %12s %12s" %('idx', 'atoms', 'ref', 'ff', 'abs diff', 'rel diff', 'weight')
        for i in xrange(len(self.pd.ric._stretches)):
            diff  = self.ref_val_str[i]-self.pd.ric.get_val_stretches()[i]
            ratio = abs(1-(self.ref_val_str[i]/self.pd.ric.get_val_stretches()[i]))*100
            print "%6d %20s %12.4f %12.4f %12.4f %12.4f %12.4f"  %\
                  (i, self.pd.ric._stretches[i], self.ref_val_str[i], self.pd.ric.get_val_stretches()[i],\
                    diff, ratio, self.wgt_val_str[i])

        print '-----------angles-----------'
        print "# %6s %20s %12s %12s %12s %12s %12s" %('idx', 'atoms', 'ref', 'ff', 'abs diff', 'rel diff', 'weight')
        for i in xrange(len(self.pd.ric._in_bends)):
            diff = (self.ref_val_ibe[i]-self.pd.ric.get_val_in_bends()[i])*rad2deg
            ratio = abs(1-(self.ref_val_ibe[i]-self.pd.ric.get_val_in_bends()[i]))
            print "%6d %20s %12.4f %12.4f %12.4f %12.4f %12.4f"  %\
                    (i, self.pd.ric._in_bends[i], self.ref_val_ibe[i]*rad2deg, self.pd.ric.get_val_in_bends()[i]*rad2deg,
                            diff, ratio, self.wgt_val_ibe[i])
        print '---------torsions-----------'
        print "# %6s %20s %12s %12s %12s %12s %12s" %('idx', 'atoms', 'ref', 'ff', 'abs diff', 'rel diff', 'weight')
        for i in xrange(len(self.pd.ric._torsions)):
            diff = (self.ref_val_tor[i]-self.pd.ric.get_val_torsions()[i])*rad2deg
            diff = 90*(signal.sawtooth(diff/rad2deg,0.5)+1) ###only multiplicity 1
            ratio = abs(1-(self.ref_val_tor[i]-self.pd.ric.get_val_torsions()[i]))
            print "%6d %20s %12.4f %12.4f %12.4f %12.4f %12.4f"  %\
                    (i, self.pd.ric._torsions[i], self.ref_val_tor[i]*rad2deg, self.pd.ric.get_val_torsions()[i]*rad2deg,
                            diff, ratio, self.wgt_val_tor[i])
        self.pd.set_xyz(self.ref(info = 'coord', branch = 'hessians', tag = self.tag))
        self.pd.write_tinker_xyz('ref.xyz')
        return



 
