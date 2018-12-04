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


class hess_fit:
    """
    class to compute the msd between the reference cartesian hessian and the models cartesian hessian.
    """
    
    def __init__(self, tag='primary', disp=0.001, weight = [0.5,0.5]):
        self.weight = np.array(weight)
        self.disp = disp
        self.tag = tag
        return

    def initialize(self, pd, ref):
        self.pd = pd
        self.ref = ref
        self.fact_hess = 1.0/143.88
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
#        if self.eigen == True:
#            self.hesscl = hessutils.hessutils(self.ref_coord, self.ref_hess,
#                    self.pd.get_elements())
#            self.ref_eigen, self.ref_evector = self.hesscl.calc_eigenvalues()
#            print np.shape(self.ref_evector)
#            print np.dot(self.ref_evector, self.ref_evector.T)
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
        self.calc_msd()
        self.cycle += 1
        self.msd = np.sum(self.weight*np.array([self.msd_coord,self.msd_hess]))
        return self.msd, [[self.msd, self.msd_coord, self.msd_hess]]
    
    def calc_msd(self):
        assert np.shape(self.coord) == np.shape(self.ref_coord)
        assert np.shape(self.hess) == np.shape(self.ref_hess)
        ### coordinational data ###
        delt_coord = self.coord - self.ref_coord
        self.msd_coord = (delt_coord*delt_coord).sum()/(self.pd.get_natoms()*3)
        ### hessian data ###
        delt_hess = (self.hess - self.ref_hess)*self.fact_hess
        self.msd_hess = (delt_hess*delt_hess).sum()/(9*self.pd.get_natoms()**2)
        self.fdiagnostics.write("%s %6.6f %6.6f\n" % (self.cycle, self.msd_coord, self.msd_hess))
        return
