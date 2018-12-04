#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
import string
import pickle
import datetime

from molsys import mpiobject
from molsys.util.timing import timer
import main


class optimizer(mpiobject):

    """
    Base class for all custom ff_gen optimizers. It is inherited from molsys mpiobject

    :Parameters:
        - ff: ff_gen class instance which should be optimized
        - mpi_comm      (mpi4py communicator, optional): Pass a local communicator, if None it uses the WORLD_COMM, default: None
        - out           (str/file object, optional): Redirect stdout to a file
    """

    def __init__(self, ff,mpi_comm=None,out=None):
        assert type(ff) == main.main
        assert ff.bsetup == True
        super(optimizer,self).__init__(mpi_comm, out)
        self.ff = ff
        self.timer = ff.timer
        return


class cma(optimizer):

    """
    CMA-ES optimizer class, inherited from optimizer and based on Niko Hansens CMA-ES optimizer
    and its python implementation. At least version 2.0.0 has to be installed. The cma class of the
    original package will is available under cma.es. The optimizer is  written in a MPI parrallel fashion.

    :Parameters:
        - ff: ff_gen class instance
        - sigma(float, optional): intial sigma, default to 0.2
        - initials(list, optionals): guessed parameters to start from, if none ff.initals will
            be used, default None
        - maxiter(int, optional): maximum number of iterations, default: 2000
        - cmaopts(dict, optional): futher options passed to the cma.es optimizer, default: {}
        - restart(str, optional): filename of a pickled cma.es instance
        - runtimelimit(str, optional): runtimelimit in the format "hh:mm:ss", default: None
    """

    def __init__(self, ff, sigma=0.2, initials = None, maxiter = 2000, cmaopts={}, restart = None, 
            runtimelimit = None):
        optimizer.__init__(self,ff)
        try:
            import cma
        except:
            raise ImportError("CMA not available --> Install via 'pip install cma'")
        assert string.split(cma.__version__,'.')[0] == '2', 'Outdated cma version installed'
        if initials is not None:
            assert len(initials)==len(ff.initials), 'Number of intial values array does not match number of variables.'
        else:
            initials = ff.initials
        self.dim = len(initials)
        options = {'bounds': ff.get_bounds(), 'maxiter':maxiter}
        options.update(cmaopts)
        self.popsize = 0
        if self.is_master:
            self.ff.pprint(41*'-'+'Optimization loop'+42*'-')
            if restart is not None:
                self.es = pickle.load(open(os.path.join('..',restart), 'rb'))
                # check if appropriate check if range is used
                assert len(initials) == len(self.es.ask()[0])
                self.es.opts.set({'maxiter':maxiter})
            else:
                self.es = cma.CMAEvolutionStrategy(initials,sigma,options)
            self.popsize = len(self.es.ask())
            self.logger = cma.CMADataLogger().register(self.es)
        self.popsize = self.mpi_comm.bcast(self.popsize, root = 0)
        # if pobj > 1 a new communicator has to be created only involving the
        # master nodes of every ffgen instance, using the Group feature of MPI
        self.idx = []
        if ff.pobj > 1:
            assert self.mpi_size % ff.pobj == 0
            mpi_group = self.mpi_comm.Get_group()
            opt_group = mpi_group.Incl(range(self.mpi_size)[::ff.pobj])
            self.opt_comm = self.mpi_comm.Create(opt_group)
        else:
            self.opt_comm = self.mpi_comm
        if self.in_optgroup:
            self.opt_size = self.opt_comm.Get_size()
            self.opt_rank = self.opt_comm.Get_rank()
            assert self.opt_size<=self.popsize
            assert self.popsize % self.opt_size == 0
#        self.popsize = self.mpi_comm.bcast(self.popsize, root = 0)
            self.nindperrank = self.popsize/self.opt_size
            self.idx = np.array(range(self.opt_rank*self.nindperrank, 
                (self.opt_rank+1)*self.nindperrank))
            # self.idx has to be communicated to ALL threads vial ff.local_comm
        self.idx = self.ff.local_comm.bcast(self.idx, root = 0)
        # handle runtimelimit
        if runtimelimit != None:
            self.runtimelimit = np.sum(np.array(map(int, string.split(runtimelimit, ":")))*np.array([3600,60,1]))
        else: 
            self.runtimelimit = np.inf
        return
    
    @property
    def in_optgroup(self):
        """
        The optgroup is am mpi group to which the master processes
        of every self.ff process belongs. If a process belongs to
        this group True is returned, else False.
        """
        return self.mpi_rank in range(self.mpi_size)[::self.ff.pobj]

    @property
    def is_master(self):
        """
        The mpi process with global rank 0 is always the master rank.
        This method returns True if the current process has rank 0, else
        False.
        """
        return self.mpi_rank == 0

    @property
    def in_mastergroup(self):
        """
        The master group is not an actual mpi group, it is a pool of
        mpiprocess which belong to the master self.ff process.
        This method returns True, if an process belongs to this pool
        else False is returned.
        """
        return self.mpi_rank in range(self.ff.pobj)

    def __call__(self, restart_interval = 100):
        """
        Method to use the optizer to optimize the actual FF. The optimizer
        stops when the stropping criterion is fullfilled or an file named
        STOP is found in the run directory.

        :Parameters:
            - restart_interval(int, optional): Interval specifying how
                often a CMA restart will be written, default: 100
        """
        self.timer.start("Fitting")
        self.converged = False
        self.result = []
        X = np.zeros([self.popsize,self.dim],dtype = float)
        Y = np.empty(self.popsize, dtype =float)
        counter = 1
        while not self.converged:
            # get solutions X
            if self.is_master:
                X = np.array(self.es.ask(),dtype=float)
            # communicate solutions to ALL threads
            self.mpi_comm.Bcast(X, root = 0)
            # calculate Y on ALL threads
            with self.timer('Generation'):
                Yr = np.array([self.ff(X[x]) for x in self.idx], dtype=float)
            # gather all solutions only from the threads in opt group to the master node
            self.mpi_comm.Barrier()
            if self.in_optgroup:self.opt_comm.Gatherv(Yr,Y, root=0)
            self.mpi_comm.Barrier()
            # tell the solutions
            if self.is_master:
                self.es.tell(X,Y)
                self.es.disp()
                self.logger.add()
                # test convergence
                if counter % restart_interval == 0:
                    pickle.dump(self.es,open('saved-cma-object.pkl','wb'))
                # check available runtime
                tpg = (datetime.datetime.now() - self.ff.starttime).total_seconds()/counter
                rest = self.runtimelimit - (datetime.datetime.now() - self.ff.starttime).total_seconds()
                if self.es.stop() or os.path.isfile("STOP") or rest < 4*tpg: self.converged = True
            # broadcast stop criterion to ALL threads
            counter += 1
            self.converged = self.mpi_comm.bcast(self.converged,root=0)
        # write output, restart and perform validations
        if self.is_master:
            self.result = self.es.result[0]
            pickle.dump(self.es,open('saved-cma-object.pkl','wb'))
        self.result = self.mpi_comm.bcast(self.result, root = 0)
        self.timer.stop("Fitting")
        self.ff.pprint(100*'-')
        if self.in_mastergroup:
            self.ff.finish(self.result)
#            self.ff.validate()
        return



