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
import matplotlib.pyplot as plt
import pickle



class energy_fit:
    """
    class to compute the msd between the reference forces and the model forces.
    """
    
    def __init__(self, tag, skip_bonded = False, ref_diff = 0.0, set_cell = False, verbose = False):
        self.tag = tag
        self.ref_diff = ref_diff
        self.skip_bonded = skip_bonded
        self.set_cell = set_cell
        self.verbose = verbose
        self.geomparams = None
        return

    def initialize(self, pd, ref):
        self.pd = pd
        self.ref = ref
        self.generate_reference()
        self.fact_energy = 1.0
        #self.get_weights()
        self.cycle = 0
        self.fdiagnostics = open('force_energy_%s.punch' % self.tag, 'w')
        return

    def set_geomparams(self, geomparams):
        self.geomparams = geomparams
        return

    def generate_reference(self):
        self.structures = copy.deepcopy(self.ref(info = 'structures',
            branch = 'forcematch', tag = self.tag))
        self.geomvectors = copy.deepcopy(self.ref(info = 'vectors',
            branch = 'forcematch', tag = self.tag))
        self.nstruc = np.shape(self.structures)[0]
        self.natoms = np.shape(self.structures)[1]
        self.energy = np.zeros((self.nstruc))
        energies = copy.deepcopy(self.ref(info = 'energies',
            branch = 'forcematch', tag = self.tag)[:,0])
        self.indmin = np.argsort(energies)[0]
        self.min_e = energies[self.indmin]
#        self.ref_energy = energies - self.min_e
        self.ref_energy = energies - self.ref_diff
        if self.set_cell == True:
            self.cells = - copy.deepcopy(self.ref(info = 'cells',
                branch = 'forcematch', tag = self.tag))
        return


    def __call__(self):
        for i in range(self.nstruc):
            if self.set_cell == True:
                self.pd.set_cell(self.cells[i,:,:])
            if self.geomparams != None:
                structure = self.structures[i,:,:] + (self.geomparams * \
                        self.geomvectors[i,:,:])
                self.pd.set_xyz(structure)
            else:
                self.pd.set_xyz(self.structures[i,:,:])
            self.pd.set_atoms_moved()
            if self.skip_bonded:
                self.pd.skip_bonded()
            e, f = self.pd.calc_energy_force()
            self.energy[i] = e
        #self.energy = self.energy-self.energy[self.indmin]
        self.calc_msd()
        self.cycle += 1
        if self.verbose:
            self.print_energies()
        return self.msd, [self.a_msd]
    
    def calc_msd(self):
        assert np.shape(self.energy) == np.shape(self.ref_energy)
        delt = (self.energy - self.ref_energy)*self.fact_energy
        self.msd = (delt*delt).sum()/self.nstruc
        self.a_msd = [self.msd]
        self.fdiagnostics.write("%s %6.6f\n" % (self.cycle, self.msd))
        return

    def print_energies(self, element = 100):
        lref = []
        lfit = []
        for i in xrange(self.nstruc):
            if i % 500 == 0:
                #lref.append(np.sum(self.ref_energy[i:i+500]+self.min_e)/self.nstruc)
                #lfit.append(np.sum(self.energy[i:i+500]+self.energy[self.indmin])/self.nstruc)
                lref.append(np.sum(self.ref_energy[i:i+500])/500.0)
                lfit.append(np.sum(self.energy[i:i+500])/500.0)
                print i, self.ref_energy[i], self.energy[i], self.ref_energy[i] -self.energy[i]
        pickle.dump(lref, open('ref.pickle', 'wb'))
        pickle.dump(lfit, open('fit.pickle', 'wb'))
        return
