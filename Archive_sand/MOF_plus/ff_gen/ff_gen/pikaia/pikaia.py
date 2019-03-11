########################################################################
##
##  pypikaia
##
##  Python wrapped pikaia derived from origianl version from
##  Marek Wojciechowski
##
##  fortran part changed to a f90 module with allocatable oldph and newph
##  arrays. generation loop moved to python with additional features to initalize
##  (read from restart, give status info etc.)
##  for use with the ff_gen force field global optimization algorithm
##
##  (C) R. Schmid, RUB, 2014
##  rochus.schmid@rub.de
##
##  derived from original Version:
##  Copyright (C) 2006 by Marek Wojciechowski
##  <mwojc@p.lodz.pl>
##
##  Distributed under the terms of the GNU General Public License (GPL)
##  http://www.gnu.org/copyleft/gpl.html
##
##
########################################################################

# import raw module
import _pikaia
# shortcut to fortran level
ppk = _pikaia.pypikaia

from mpi4py import MPI
import numpy as np
import time
import string

class pypikaia:

    """
    Pypikaia (Pythonified F90 Module version of 
    Pikaia V1.2 (genetic algorithm based optimizer)
    
    Main features:
        - Access to values from python
        - Writes and reads restart files
        - Allows to set variable ranges for each parameter
        - Generation loop on python level

    :Parameters:
        - func : callable
            Scalar function of the signature fitness = func(x),
            where *x* is a real array of length *n*.
            By default *x* must be in the range [0,1] but with 
            the method set_bounds other ranges can be defined
            By convention, ff should return higher fitness values for more
            optimal parameter values (i.e., individuals which are
            more "fit"). For example, in fitting a function
            through data points, *func* could return the inverse of
            chi**2.
        - n : int
            Length of *x*. 
        - individuals : int
            Number of individuals in a population (default is 100)
        - digits :  int
            Number of significant digits (i.e., number of
            genes) retained in chromosomal encoding (default
            is 6).
        - crossover :  float
            Crossover probability; must be  <= 1.0 (default
            is 0.85). If crossover takes place, either one
            or two splicing points are used, with equal
            probabilities.
        - mutation : {1, 2, 3, 4, 5, 6}
            =====   =====================================================
            digit   description
            =====   =====================================================
              1     one-point mutation, fixed rate
              2     one-point, adjustable rate based on fitness (default)
              3     one-point, adjustable rate based on distance
              4     one-point+creep, fixed rate
              5     one-point+creep, adjustable rate based on fitness
              6     one-point+creep, adjustable rate based on distance
            =====   =====================================================
        - initrate : float
            Initial mutation rate. Should be small (default
            is 0.005) (Note: the mutation rate is the probability
            that any one gene locus will mutate in any one generation.)
        - minrate : float
            Minimum mutation rate. Must be >= 0.0 (default is 0.0005)
        - maxrate : float
            Maximum mutation rate. Must be <= 1.0 (default is 0.25)
        - fitnessdiff : float
            Relative fitness differential. Range from 0 (none)
            to 1 (maximum) (default is 1).
        - reproduction : {1, 2, 3}
            Reproduction plan; 1/2/3 = Full generationalreplacement/Steady-
            state-replace-random/Steady-state-replace-worst (default is 3)
        - elitism : {0, 1}
            Elitism flag; 0/1=off/on (default is 0)
            (Applies only to reproduction plans 1 and 2)
        - verbosity : {0, 1, 2}
            Printed output 0/1/2=None/Minimal/Verbose (default is 0)
        - true_random : {True / False} (default True)
            Initalize the random generator to a give different numbers each run
            Note that the seed is broadcasted from the Master node so that this can be used
            also in MPI parallel runs (each node will produce the same random numbers)

    .. note::
        Original fortran code of pikaia is written by:
        Paul Charbonneau & Barry Knapp (paulchar@hao.ucar.edu,
        knapp@hao.ucar.edu)

        Original Python wrapped version by Marek Wojciechowski (mwojc@p.lodz.pl)
        
        Extended F90-Module based version with PY-level generation loop by Rochus Schmid (rochus.schmid@rub.de)
    """

    
    def __init__(self, func, n, \
            individuals = 100, \
            digits = 6, \
            crossover = 0.85, \
            mutation = 2, \
            initrate = 0.005, \
            minrate = 0.0005, \
            maxrate = 0.25, \
            fitnessdiff = 1.0, \
            reproduction = 3, \
            elitism = 0, \
            verbosity = 0, \
            true_random=True):

        self.comm = MPI.COMM_WORLD 
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        # set control variables
        self.n  = n
        self.func = func
        self.np = individuals
        self.nd = digits
        # set fortran level control variables
        ppk.pcross = crossover
        ppk.imut   = mutation
        ppk.pmut   = initrate
        ppk.pmutmn = minrate
        ppk.pmutmx = maxrate
        ppk.fdif   = fitnessdiff
        ppk.irep   = reproduction
        ppk.ielite = elitism
        ppk.ivrb   =verbosity
        # if true random is used initalize the random generator from numpy
        if true_random:
            seed = np.random.randint(999999)
            if self.size > 1:
                seed = self.comm.bcast(seed, root=0)
        else:
            seed = 0
        # initalize and allocate arrays and check control variables
        status = ppk.setup(self.n, self.np, self.nd, seed)
        assert status == 0, "Error during setup!!!"
        # set some defaults
        self.minbound = np.zeros([self.n],dtype="float64")
        self.maxbound = np.ones([self.n],dtype="float64")
        self.use_bound=False
        self.ngen = 0
        return
        
    def set_bounds(self, minbound, maxbound):
        """
        set boundary conditions for each parameter
        
        Internally pikaia treats parameters between 0 and 1
        If this function is called the values are used to convert the paramters from the 
        given range into values between 0 an 1
        
        :Paramter:
             - minbound : list or numpy array len(n) with lower bound
             - maxbound : list or numpy array len(n) with upper bound
        """     
        
        assert len(minbound) == self.n
        assert len(maxbound) == self.n
        self.minbound = np.array(minbound, dtype="float64")
        self.maxbound = np.array(maxbound, dtype="float64")
        assert (self.minbound < self.maxbound).all()
        self.valrange = self.maxbound-self.minbound
        self.use_bound = True
        return
        
    def convert_from_internal(self, params):
        """
        converts from internal range to specifeid range if bounds are used
        """
        if self.use_bound:
            params = self.minbound + (self.valrange*params)
        return params
        
    def convert_to_internal(self, params):
        """
        converts from internal range to specifeid range if bounds are used
        """
        if self.use_bound:
            params = (params - self.minbound)/self.valrange
        return params
        
    def invalidate_fitness(self):
        """
        calling this means the fittness of the population needs to be computed and ranked again
        """
        ppk.ftns_calculated = False
        ppk.ftns_sorted     = False
        return
    
    def invalidate_ranking(self):
        """
        calling this means the ranking of the fittness needs is not valid anymore
        (for example if population and fittness have been read from file)
        """
        ppk.ftn_sorted      = False
        return
        
    def get_all_values(self):
        """
        get values from fortran array
        in sorted form
        
        :Returns:
            - vars : list of list with parameter values (sorted by fittness)
            - fits : list of fittness
        """
        ilist = ppk.ifit.tolist()
        ilist.reverse()
        fitlist = []
        vallist = []
        for i in ilist:
            val = ppk.oldph[:,i-1]
            fit = ppk.fitns[i-1]
            val = self.convert_from_internal(val)
            fitlist.append(fit)
            vallist.append(val)
        return vallist, fitlist

    def set_all_values(self, values, fitns):
        """
        set all values and fitns to some intial values
        Ranking is invalidated
        
        :Parameters:
            - values : numpy array of shape [population, size]
            - fitns  : numpy array of shape [population]
        """
        values = self.convert_to_internal(values)
        ppk.oldph[:,:] = values.transpose()
        ppk.fitns[:]   = fitns
        self.invalidate_ranking()
        return

    def set_value(self, values, fitns, ind):
        """
        set values and fitns for a single species to some intial values
        Ranking is invalidated
        
        :Parameters:
            - values : numpy array of shape [size]
            - fitns  : float
            - ind    : integer, index of species (between 0 and population)
        """
        assert (ind >= 0) and (ind < self.np)
        values = self.convert_to_internal(values)
        ppk.oldph[:,ind] = values
        ppk.fitns[ind]   = fitns
        self.invalidate_ranking()
        return

    def get_all_fitness(self):
        """
        get fitness in sorted form and average fittness
        
        :Returns:
            - fits : list of fittness in sorted form
            - aver : average value
        """
        fitlist = np.sort(ppk.fitns)
        aver = np.average(fitlist)
        fitlist = fitlist.tolist()
        fitlist.reverse()
        return fitlist, aver

    def callback(self, params):
        """
        callback passed to fortran
        converts parmaeters if necessary .. do not call this method
        """
        params = self.convert_from_internal(params)
        fitness = self.func(params)
        return fitness
                
    def initialize(self):
        """
        initalize fitness and rank it if necessary
        :Note:
        This is controlled by the flags ftns_calculated and ftns_sorted.
        If fitness is once computed and ranked these stay true because any GA cycle
        automatically computes fittnes and ranks it for a new population.
        So this needs only to be done at startup or if any other manipulation of the population
        is performed.
        Use the methods invalidate_fitness and invalidate_ranking to trigger this
        """
        ppk.initialize(self.callback, self.n, self.np)
        return
        
    def iterate(self, steps, thresh = None, restart=None):
        """
        Do a number of GA iteration starting from the current poplation
        (initial fittness and ranking is performed if necessary)
        
        :Parameters:
            - steps   : Number of steps to run
            - restart : Name of the restart file or *None* (default, no restart is written)
            - thresh  : threshold of fittness to terminate (if this fittness is reached by the best stop)
        """
            
        self.initialize()
        if self.rank == 0 : 
            print ("Performing %5d GA iterations" % steps)
            ff = open('pikaia.out', 'w')
            fit, aver = self.get_all_fitness()
            ff.write("0 best: %10.5f  worst: %10.5f  average: %10.5f\n" % (fit[0], fit[-1], aver))
            print ("0 best: %10.5f  worst: %10.5f  average: %10.5f" % (fit[0], fit[-1], aver))
        for i in xrange(steps):
            self.ngen += 1
            status = ppk.pikaia(self.callback, self.n, self.np, self.nd, self.ngen)
            fit, aver = self.get_all_fitness()
            #
            #values = self.get_all_values()
            if self.rank == 0 : 
                print "Generation %d" % (self.ngen)
                print ("%d best: %10.5f  worst: %10.5f  average: %10.5f" % (self.ngen, fit[0], fit[-1], aver))
                ff.write("%d best: %10.5f  worst: %10.5f  average: %10.5f\n" % (self.ngen, fit[0], fit[-1], aver))
            if restart: self.write_restart(restart)
            if thresh:
                if fit[0] >= thresh:
                    if self.rank == 0:
                        print "Requested fittness has been reached! Stopping iterations."
                    break
        if self.rank == 0:
            ff.close()
        return
        
    def write_restart(self, filename):
        """
        wrtie restart information to file
        
        :Paramters:
            - filename : string
        """
        if self.rank == 0:
            rf = open(filename, "w")
            rf.write("# pypikaia restart file V1.0\n")
            rf.write("# %s\n" % time.asctime())
            rf.write("Size: %5d  Population: %5d  Generation: %5d\n" % (self.n, self.np, self.ngen))
            val, fit = self.get_all_values()
            for i,v in enumerate(val):
                rf.write("%5d : %15.8f \n" % (i+1, fit[i]))
                for j in xrange(self.n):
                    rf.write("%15.8f " % v[j])
                rf.write("\n")
            rf.close()
        return
    
    def read_restart(self, filename):
        """
        read a complete population from a restart file 
        
        :Parameters:
            - filename : filename of the restart file (V1.0)
        """
        rf = open(filename, "r")
        l = string.split(rf.readline())
        assert l[4] == "V1.0" , "Wrong version of restart file"
        l = string.split(rf.readline())
        if self.rank == 0 : print("Reading restart file %s from %s" % (filename,string.join(l[1:])))
        l = string.split(rf.readline())
        n   = int(l[1])
        npop  = int(l[3])
        ngen  = int(l[5])
        assert self.n  == n
        assert self.np == npop
        self.ngen = ngen
        fit = np.zeros([self.np], dtype="float64")
        val = np.zeros([self.np, self.n], dtype="float64")
        for i in xrange(self.np):
            l = string.split(rf.readline())
            assert i+1 == int(l[0])
            fit[i] = float(l[2])
            l = string.split(rf.readline())
            v = np.array(map(string.atof, l))
            val[i,:] = v
        rf.close()
        self.set_all_values(val, fit)
        return
            
    def read_multiple_restart(self, filelist, nindividuals=None):
        """
        read more than one restart file given in a list
        From each file are the best *nindividuals* read in.
        If this is not specified then nindividuals is determined from the size of then
        population and the number of files given trying to fill the polpulation complete
        Any remaining species are initalized randomly
        
        :Parameters:
            - filelist : List of files from which restart species will be read in
            - nindividuals : number of individuals to be read from each file 
            
        """
        nfiles = len(filelist)
        if not nindividuals:
            nindividuals = self.np/nfiles
        assert nindividuals*nfiles <= self.np, "Population is too small to read so many individuals"
        ind = 0
        for f in filelist:
            rf = open(f, "r")
            l = string.split(rf.readline())
            assert l[4] == "V1.0" , "Wrong version of restart file"
            l = string.split(rf.readline())
            if self.rank == 0 : print("Reading restart file %s from %s" % (restart,string.join(l[1:])))
            l = string.split(rf.readline())
            n   = int(l[1])
            npop  = int(l[3])
            ngen  = int(l[5])
            assert self.n  == n
            assert npop > nindividuals
            for i in xrange(nindividuals):
                l = string.split(rf.readline())
                fit = float(l[2])
                l = string.split(rf.readline())
                v = np.array(map(string.atof, l))
                self.set_value(v, fit, ind)
                ind += 1
            rf.close()
        if ind < self.np:
            if self.rank == 0: print("A total of %d species were read from file, filling up randomly")
            nrand = self.np-ind
            for i in xrange(nrand):
                randpar = np.random.random([self.n])
                if self.size>1:
                    randpar = self.comm.bcast(randpar, root=0)
                # since randpar is in the 0-1 range we can directly set to oldph insread of using set_values
                ppk.oldph[:,ind] = randpar
                # now compute fittness directly here .. 
                fitns = self.callback(randpar)
                ppk.fitns[ind] = fitns
                ind += 1
        self.ngen = 0
        self.invalidate_ranking()
        return
                
    
    def random_start(self):
        """
        start up with random parameters
        """
        randpar = np.random.random([self.n, self.np])
        if self.size > 1:
            randpar = self.comm.bcast(randpar, root=0)
        ppk.oldph[::] = randpar
        self.invalidate_fitness()
        return

    
