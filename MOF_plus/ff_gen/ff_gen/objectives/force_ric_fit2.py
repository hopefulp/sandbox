"""
This file implements a ric_fit class.
It is inheriting all features from RedIntCoords in ric.py
and adds features to load a respective reference structure and Hessian.
In addition, a weight matrix is held to evaluate various kinds of weighted mean
square deviations to be used as ingredients to fittness values
"""

import string
import numpy as np
import copy
import os

from ff_gen.objectives.base import base
#from ff_gen.refclass import refclass2 as refclass
from ff_gen.objectives.ric_fit3 import mapper

#import pydlpoly
#import pylmps



class force_ric_fit(base):
    """
    class to compute redundant internal coordinates (ric)
    by using the inherited ric module and to compute deviations from a corresponding
    reference.
    """
    
   
    def __init__(self, pd, reffile, tag, start_dir, fullric = True, fact_force = 1.0,mpi_comm=None,out=None):
        super(force_ric_fit,self).__init__(pd,reffile,tag,start_dir,mpi_comm, out)
        self.name = "ForceRicFit"
        self.fact_force = fact_force
        self.pd.mol.addon('ric')
        self.pd.mol.ric.setup_rics(full=fullric)
        self.generate_reference(self.pd.get_cell())
        self.set_weights()
        #self.initialize(fullric)
        return

    def set_weights(self, norm = True):
        """ Set the weights for the individual rics """
        varpots    = self.pd.mol.ff.par.variables.varpots
        varpotnums = self.pd.mol.ff.par.variables.varpotnums(self.pd.mol.ff)
        ics = ["bnd", "ang", "dih", "oop"]
        self.wgt = np.zeros([self.nstruc, self.pd.mol.ric.num_ric])
        for ic in ics:
            for i,r in enumerate(self.pd.mol.ff.ric_type[ic]):
                j, j_glob = self.pd.mol.ric.map_ric(mapper[ic],r)
                ### check for linear bend case
                if ic == "ang" and r.lin == True: 
                    j, j_glob = self.pd.mol.ric.map_ric("lbe",r)
                if j == None: continue 
                pi = self.pd.mol.ff.parind[ic][i]
                weight = 0.0
                for p in pi:
                    if (ic, p) in varpots:
                        weight = 1.0
                        if norm: weight/=varpotnums[(ic, p)]
                        break
                if ic == "oop":
                    self.wgt[:,j_glob]   = weight/3.0
                    self.wgt[:,j_glob+1] = weight/3.0
                    self.wgt[:,j_glob+2] = weight/3.0
                    #import pdb; pdb.set_trace()
                elif ic == "ang" and r.lin == True:
                    self.wgt[:,j_glob]   = weight/2.0
                    self.wgt[:,j_glob+1] = weight/2.0
                else:
                    self.wgt[:,j_glob] = weight
        return

    def project_force(self, cforce, invb):
        """ Project force from Cart to RICs """
        rforce = np.dot(np.transpose(invb), np.ravel(cforce))
        return rforce

    def generate_reference(self, cell):
        """
        Method to generate the reference data
        """
        # get structures out of the reference set
        self.structures = copy.deepcopy(self.ref(info = 'coords', branch = 'forcematch', tag = self.tag))
        # get nstruc
        self.nstruc = np.shape(self.structures)[0]
        # get cartesian refererence forces
        self.ref_cart_force = -1.0*copy.deepcopy(self.ref(info = 'gradients', branch = 'forcematch', tag = self.tag))
        # allocate ref_force arrays storing the force in rics
        self.ref_ric_force = np.zeros([self.nstruc, self.pd.mol.ric.num_ric])
        self.ric_force = np.zeros([self.nstruc, self.pd.mol.ric.num_ric])
        self.cart_force = np.zeros([self.nstruc, self.pd.mol.natoms, 3])
        # list of inv_bmats needed to project the forces into rics
        self.inv_bmats = []
        for i in range(self.nstruc):
            # construct b matrix
            self.pd.mol.ric.construct_b_matrix(cell, self.structures[i,:,:])
            # ivert b matrix
            invb, rankb = self.pd.mol.ric.invert_b_matrix()
            if rankb < self.pd.mol.natoms: raise ValueError('Rank too small!')
            self.inv_bmats.append(invb)
            # project force
            self.ref_ric_force[i,:] = copy.deepcopy(self.project_force(self.ref_cart_force[i,:,:],
                invb))
        return
    
    def __call__(self):
        """
        Returns the fitness(msd) of the current MM model
        """
        # write variables to the calculators internal data structure
        self.writer(self.pd)
        # loop over all structures and calculate the MM force
        # project it and then calc the msd to the reference
        for i in range(self.nstruc):
            self.pd.set_xyz(self.structures[i,:,:])
            self.pd.set_atoms_moved()
            e, f = self.pd.calc_energy_force()
            self.ric_force[i,:] = self.project_force(f, self.inv_bmats[i])
            self.cart_force[i,:,:] = f
        self.cycle += 1
        return self.msd

    @property
    def msd(self):
        """
        Calculates the mean square deviations and returns them as a np array
        """
        wdelt = (self.ric_force-self.ref_ric_force)*self.fact_force*self.wgt
        msd = (wdelt*wdelt).sum()/self.nstruc
        return msd
     
    def finish(self):
        """ 
        Method is called after the end of the fitting process to print analysis
        and performance files
        """
        retval = os.getcwd()
        os.chdir(self.pd.rundir)
        ### your analysis methods are called here 
        os.chdir(retval)

