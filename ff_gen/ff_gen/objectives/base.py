#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import os
import sys

from mpi4py import MPI

from molsys import mpiobject
from ff_gen.refclass import refclass2 as refclass

class base(mpiobject):

    """
    Base class ob every objective class which should be used in the FFgen
    framework, class is inherited from molsys.mpiobject

    :Parameters:
        - pd  (obj): initialized pydlpoly or lammps object used as MM machine
        - ref (str): name of the reffile
        - tag (str): Name of the subgrpup in the reference file holding the
            reference information
        - start_dir(str): start_dir of the ffgen launch
        - mpi_comm(obj, optional): MPI communicator which should be used, defaults to None
        - out(str, optional): Name of the out file, which should be used for printing, defaults to None
    """


    def __init__(self, pd, reffile, tag, start_dir, mpi_comm = None, out = None):
        super(base, self).__init__(mpi_comm, out)
        assert type(tag)  == str
        assert type(reffile)  == str
        assert type(start_dir) == str
        self.name = 'Base'
        self.pd  = pd
#        import pdb; pdb.set_trace()
        if pd.__class__.__name__ == "pylmps":
            self.backend = 'lammps'
            self.writer = self.pd.ff2lmp.write2internal
        elif pd.__class__.__name__ == "pydlpoly":
            self.backend = 'pydlpoly'
            self.writer = self.pd.mol.write2internal
        else:
            raise NotImplementedError('Unkown MM backend passed!')
        self.ref = refclass(reffile)
        self.tag = tag
        self.start_dir = start_dir
        self.cycle = 0
        return

    def __repr__(self):
        return "%s-%s-%s" % (self.name, self.ref.name, self.tag)


    def __call__(self):
        """
        Method called by the ff_gen class to evaluate the fitness
        of a parameter set, should return the fitness.
        """
        raise NotImplementedError('This function has to be overloaded')


    def finish(self):
        """
        Method called after finishing a parameterization run by the ff_gen class.
        All analysis after a run has finished should go here.
        """
        ### call after the last __call__ with final params
        ### move to rundir of the calculator
        if self.mpi_rank != 0: return
        retval = os.getcwd()
        os.chdir(self.pd.rundir)
        # write optimized params
        self.pd.mol.ff.write('opt')
        return

