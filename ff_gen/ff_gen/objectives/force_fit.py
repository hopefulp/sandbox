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


class force_fit:
    """
    class to compute the msd between the reference forces and the model forces.
    """
    
    def __init__(self, tag, cenergy = False, weights = None, skip_bonded = False, set_cell = False):
        self.tag = tag
        self.cenergy = cenergy
        self.skip_bonded = skip_bonded
        self.set_cell = set_cell
        if weights == None:
            self.weights = [0.5,0.5]
        else:
            self.weights = weights
        return

    def initialize(self, pd, ref):
        self.pd = pd
        self.ref = ref
        self.fact_force = 1.0/143.88
        self.generate_reference()
        self.get_weights()
        self.cycle = 0
        self.fdiagnostics = open('force_fit_%s.punch' % self.tag, 'w')
        return

    def generate_reference(self):
        self.structures = copy.deepcopy(self.ref(info = 'structures',
            branch = 'forcematch', tag = self.tag))
        self.nstruc = np.shape(self.structures)[0]
        self.natoms = np.shape(self.structures)[1]
        self.force = np.zeros([self.nstruc, self.natoms,3], dtype = 'float64')
        self.energy = np.zeros([self.nstruc], dtype = 'float64')
        self.wgt_force = np.zeros([self.nstruc, self.natoms,3], dtype = 'float64')
        self.ref_force = - copy.deepcopy(self.ref(info = 'forces',
            branch = 'forcematch', tag = self.tag))
        if self.cenergy == True:
            energies = copy.deepcopy(self.ref(info = 'energies',
                branch = 'forcematch', tag = self.tag)[:,0])
            min_e = energies[0]
            self.index_min_e = 0
            for i in range(self.nstruc - 1):
                if energies[i+1] < min_e:
                    min_e = energies[i+1]
                    self.index_min_e = i+1
            self.ref_energy = energies - min_e
        if self.set_cell == True:
            self.cells = - copy.deepcopy(self.ref(info = 'cells',
                branch = 'forcematch', tag = self.tag))
#        print 'force:', self.ref_force[:,0,0]
#        print 'force:', self.ref_force[:,1,1]
#        print 'force:', self.ref_force[:,2,2]
        print np.sum(np.abs(self.ref_force), axis = (1,2))
#        print self.ref_force
#        print self.wgt_force
#        print self.cells
#        print self.ref_force
        return

    def get_weights(self):
#        self.pd.mol.var_atoms = [0,1]
        for i in range(self.natoms):
            if i in self.pd.mol.var_atoms:
                self.wgt_force[:,i,:] = 1.0
        return

    def __call__(self):
        for i in range(self.nstruc):
            if self.set_cell == True:
                self.pd.set_cell(self.cells[i,:,:])
            self.pd.set_xyz(self.structures[i,:,:])
            self.pd.set_atoms_moved()
            if self.skip_bonded:
                self.pd.skip_bonded()
            e, f = self.pd.calc_energy_force()
            self.force[i,:,:] = f
            self.energy[i] = e
        if self.cenergy:
            self.energy = self.energy-self.energy[self.index_min_e]
        self.calc_msd()
        self.cycle += 1
        return self.msd, [self.a_msd]
    
    def calc_msd(self):
        assert np.shape(self.force) == np.shape(self.ref_force)
        delt = (self.force - self.ref_force)*self.fact_force
        wdelt = delt * self.wgt_force
        self.f_msd = (wdelt*wdelt).sum()/self.wgt_force.sum()
        if self.cenergy == True:
            edelt = (self.energy - self.ref_energy)
            self.e_msd  = (edelt*edelt).sum()/self.nstruc
            self.msd = self.weights[0] * self.f_msd +\
                    self.weights[1]*self.e_msd
            self.a_msd = [self.msd, self.f_msd, self.e_msd]
            self.fdiagnostics.write("%s %6.6f %6.6f\n" % (self.cycle, self.f_msd, self.e_msd))
        else:
            self.msd = self.f_msd
            self.a_msd = [self.msd]
            self.fdiagnostics.write("%s %6.6f\n" % (self.cycle, self.msd))
        return
