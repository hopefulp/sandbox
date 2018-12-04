"""
This file implements a force_fit class.
It handles the reference force arrays and calculates the msd.
In addition, a weight matrix is held to evaluate various kinds of weighted mean
square deviations to be used as ingredients to fittness values
"""

import string
import numpy as np
import copy
import os
import refclass
import hessutils


class freq_fit:
    """
    class to compute the msd between the reference nm freqs and the models nm freqs.
    """
    
    def __init__(self, tag='primary', disp=0.001, weight = [0.5,0.5], freq_weight = None):
        self.weight = np.array(weight)
        self.disp = disp
        self.tag = tag
        self.freq_weight = freq_weight
        return

    def initialize(self, pd, ref):
        self.pd = pd
        self.ref = ref
        self.fact_freq = 1.0
        if self.freq_weight == None:
            self.freq_weight = np.ones((self.pd.get_natoms()*3))
        else:
            assert len(self.freq_weight) == self.pd.get_natoms()*3
        self.generate_reference()
        #self.get_weights()
        self.cycle = 0
        self.fdiagnostics = open('hess_fit_%s.punch' % self.tag, 'w')
        return

    def generate_reference(self):
        self.ref_coord = copy.deepcopy(self.ref(info='coord', branch='hessians',
            tag = self.tag))
        self.ref_hess = copy.deepcopy(self.ref(info='hessian', branch='hessians',
            tag = self.tag))
        self.hesscl = hessutils.hessutils(self.ref_coord, self.ref_hess,
                self.pd.get_elements())
        self.ref_eigen, self.ref_evector = self.hesscl.calc_eigenvalues()
        print 'freqs'
        print self.ref_eigen
        return

#    def get_weights(self):
#        self.pd.mol.var_atoms = [0,1]
#        for i in range(self.natoms):
#            if i in self.pd.mol.var_atoms:
#                self.wgt_force[:,i,:] = 1.0
#        return

    def __call__(self):
        self.pd.MIN_lbfgs(0.001)
        self.coord = self.pd.get_xyz()
        self.hess = self.pd.calc_hessian(disp = self.disp)
        self.hesscl.s_hessian = self.hess
        self.hesscl.xyz = self.coord
        self.eigen, self.evector = self.hesscl.calc_eigenvalues()
        self.calc_msd()
        self.cycle += 1
        self.msd = np.sum(self.weight*np.array([self.msd_coord,self.msd_freq]))
        return self.msd, [[self.msd, self.msd_coord, self.msd_freq]]
    
    def calc_msd(self):
        assert np.shape(self.coord) == np.shape(self.ref_coord)
        assert np.shape(self.hess) == np.shape(self.ref_hess)
        ### coordinational data ###
        delt_coord = self.coord - self.ref_coord
        self.msd_coord = (delt_coord*delt_coord).sum()/(self.pd.get_natoms()*3)
        ### frequency data ###
        delt_eigen = (self.eigen - self.ref_eigen)*self.freq_weight*self.fact_freq
        self.msd_freq = (delt_eigen**2).sum()/(np.sum(self.freq_weight)*3)
        self.fdiagnostics.write("%s %6.6f %6.6f\n" % (self.cycle, self.msd_coord, self.msd_freq))
        return
