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
import ff_gen.refclass as refclass
import matplotlib.pyplot as plt
import pickle



class grid_fit:
    """
    class to compute the msd between the reference forces and the model forces.
    """
    
    def __init__(self, tag, skip_bonded = False, ref_diff = 0.0, set_cell = False, verbose = False, iatom = -1):
        self.tag = tag
        self.ref_diff = ref_diff
        self.skip_bonded = skip_bonded
        self.set_cell = set_cell
        self.verbose = verbose
        self.geomparams = None
        self.iatom = iatom
        return

    def initialize(self, pd, ref, gen_ref = True):
        self.pd = pd
        self.ref = ref
        if gen_ref == True:
            self.generate_reference()
        self.fact_energy = 1.0
        self.cycle = 0
        #self.fdiagnostics = open('force_energy_%s.punch' % self.tag, 'w')
        return

    def eval_grid(self, grid, iatom):
        energies = []
        for i in range(np.shape(grid)[0]):
            xyz = self.pd.get_xyz()
            xyz[iatom,:] = grid[i,:]
            self.pd.set_xyz(xyz)
            self.pd.set_atoms_moved()
            e, f = self.pd.calc_energy_force()
            energies.append(e)
        energies = np.array(energies)
        return energies

    def generate_reference(self):
        self.grid = copy.deepcopy(self.ref(info = 'grid', 
            branch = 'grids', tag = self.tag))
        self.ngrid = copy.deepcopy(self.ref(info = 'ngrid',
            branch = 'grids', tag = self.tag))
        self.energies = copy.deepcopy(self.ref(info = 'energy',
            branch = 'grids', tag = self.tag))
        return


    def __call__(self):
        energies = self.eval_grid(self.grid, self.iatom)
        self.msd = np.sum((energies - self.energies)**2)/self.ngrid[0,0]
        self.cycle += 1
        return self.msd, [self.msd]
