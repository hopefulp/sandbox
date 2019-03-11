# -*- coding: utf-8 -*- 
import os
import atexit
import getpass
import socket
import copy
import shutil
import string
import datetime
import numpy as np
import mpi4py
import subprocess

# try to import the mm machines, catch if only one is installed
mmclasses = {}
try:
    import pydlpoly
    mmclasses['pydlpoly']=pydlpoly.pydlpoly
except ImportError:
    pass
try:
    import pylmps
    mmclasses['lammps']=pylmps.pylmps
except ImportError:
    pass
assert len(mmclasses) > 0, "No MM backend available! Please install PYDLPOLY or PYLMPS."


import molsys
from molsys.util.ffparameter import potentials
from molsys.util.timing import timer, Timer

import ff_gen
from objectives import ric_fit3 as ric_fit
from objectives import force_ric_fit2 as force_ric_fit
#pydlpoly.dlp_excl.mapper[3] = False

# setup logger object
import logging
logger = logging.getLogger("FFgen")
logger.setLevel(logging.INFO)
fhandler = logging.FileHandler("FFgen.log")
fhandler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m-%d %H:%M')
fhandler.setFormatter(formatter)
logger.addHandler(fhandler)

def find_rundir(name, rundir = None):
    """
    Method to find a rundir for a given name, if directrory already exist it will increment the name.

    :Parameters:
        - name(str): name of directory
        - rundir(str, optional): if not None use this directory, if it already exist raise IOError,
            default: None

    :Returns:
        - temprundir(str): final directory name
    """
    if rundir is not None:
        if os.path.isdir(rundir): 
            raise IOError('rundir %s already there' % rundir)
        else:
            return rundir
    else:
        temprundir = name
        i = 1
        while os.path.isdir(temprundir):
            i+=1
            temprundir = '%s_%d' % (name, i)
        return temprundir


class main(molsys.mpiobject):
    """
    The ff_gen class is the main class of the FFgen code. It makes the actual
    objective function available to arbitrary optimizer. It is inherited from 
    molsys.mpiobject. If the mpi_size is larger than one the communicator is
    splitted in a way that every global rank gets its own local comm.

    :Parameters:
        - name          (str) : name of the ff_gen instance, a directory with this name will be created if rundir is None 
        - minimize      (bool, optional): flag to toggle from maximizing to minimizing the objective function, default: False
        - use_range     (bool, optional): if used the parameters are scaled between 0 and 1 by using their indivual ranges, default: False
        - rundir        (str, optional): by default the rundir is created based on the name parameter, if an explicit name should be used
            rundir has to be used, default: None
        - pobj          (int, optional): Integer specifying into how many local communicators the global MPI communicator is splitted, default: 1
        - mpi_comm      (mpi4py communicator, optional): Pass a local communicator, if None it uses the WORLD_COMM, default: None 
        - out           (str/file object, optional): Redirect stdout to a file
    """

    def __init__(self, name, minimize = True, use_range = False, rundir = None,
            pobj = 1, mpi_comm = None, out = None):
        super(main,self).__init__(mpi_comm, out)
        # set name
        self.name = name
        # do mpi stuff
        assert type(pobj) == int, "pobj has to be an integer"
        assert pobj >= 1, "pobj hast to be at least 1"
        assert self.mpi_size % pobj == 0, "not possible to distribute %i parallel objectives over %i cores" % (pobj, self.mpi_size)
        ### check mpi4py version
        if self.mpi_size > 1 and self.mpi_rank == 0:
            if mpi4py.__version__ != "2.0.0": 
                raise NotImplementedError("MPI parrallel FFgen is only working with mpi4py version 2.0.0")
        self.pobj = pobj
        if self.mpi_size > 1:
            if pobj == 1:
                self.local_comm = self.mpi_comm.Split(self.mpi_rank, 0)
            else:
                self.local_comm = self.mpi_comm.Split(self.mpi_rank/2,self.mpi_rank)
        else:
            self.local_comm = self.mpi_comm
        self.local_rank = self.local_comm.Get_rank()
        self.local_size = self.local_comm.Get_size()
        # create directory and move to it
        self.start_dir = os.getcwd()
        self.rundir = find_rundir(self.name, rundir)
        self.mpi_comm.Barrier()
        if self.mpi_rank == 0: os.mkdir(self.rundir)
        self.mpi_comm.Barrier()
        os.chdir(self.rundir)
        self.print_startup()
        # init parameter dictionary
        self.par = potentials({
                    "cha": {},
                    "vdw": {},
                    "bnd": {},
                    "ang": {},
                    "dih": {},
                    "oop": {},
                    })
        self.par.FF=name
        self.par.attach_variables()
        # further settings
        self.minimize = minimize
        self.use_range = use_range
        self.cycle = 0
        self.objs  = []
        self.weights = []
        self.bsetup = False
        self.finished = False
        self.validations = []
        atexit.register(self.print_finish)
        self.timer = Timer()
        return

    def status_print(self, stype, statement):
        """
        Helper function to make output status prints in the CP2K style
        
        :Parameters:
            - stype (str): General keyword
            - statement (str): Actual statement
        """
        self.pprint("%s|%s" % (stype ,statement))
        return

    def print_startup(self):
        """
        Method called at startup to print general info.
        """
        if self.mpi_rank != 0: return
        self.starttime = datetime.datetime.now()
#        rev,branch = molsys.util.hg.get_revison()
        self.pprint(ff_gen.header)
        self.status_print('FFgen', "version string: %s" % ff_gen.__version__)
        self.status_print('Global',"Total number of message passing processes: %i" % self.mpi_size)
        self.status_print('Global',"Hostname: %s" % socket.gethostname())
        self.status_print('Global',"User: %s" % getpass.getuser())
        self.status_print('Global',"Starttime: %s" % self.starttime)
        self.status_print('Global',"Startdir: %s" % self.start_dir)
        self.status_print('Global',"Rundir: %s" % os.path.join(self.start_dir,self.rundir))
        return

    def print_finish(self):
        """
        Method called via atextit to print general info incl timings.
        """
        if self.mpi_rank != 0: return
        self.status_print('Global', 'Endtime: %s' % datetime.datetime.now())
        self.timer.write(date=False)
#        self.timer.write_logger(logger.info)
        self.pprint(ff_gen.footer)
        # move back to start dir
        os.chdir(self.start_dir)
        return

    @timer("init_objective")
    def init_objective(self, objname, sysname, tag, coord=None, fpar=None, reffile=None, weight = 1.0, offline = False, 
            refsys = None, mm = 'lammps', bcond =3, stdout = False,objargs = [], objkwargs = {}, mmkwargs = {}):
        """
        Method to initialize an objective class
        
        :Parameters:
            - objname   (str): name of the objective class
            - sysname   (str): Name of the molecular system, if coord or fpar or reffile are None, the sysname is
                also used for this files
            - tag       (str): Name of the subgroup in the hdf5 file holding the actual reference information
            - coord     (str, optional): Name of the mfpx file, if None sysname.mfpx will be used, default: None
            - fpar      (str, optional): Name of the fpar/ric file, if None sysname.fpar will be used, default: None
            - reffile   (str, optional): Name of the hdf5 reference file, if None sysname.hdf5 will be used, default: None
            - weight    (float, optional): Weight of the objective in the total objective function, default: 1.0
            - offline   (bool, optional): Flag to specify if molsys offline assignment should be used, default: False
            - refsys    (str, optional): Name of the reference system in case of offline assignment, default: None
            - mm        (str, optional): Name of the molecular mechanics engine to use, only pydlpoly or lammps possible, default: lammps
            - bcond     (int, optional): Bconds used in MM engine, default: 3
            - stdout    (bool, optional): Flag to toggle if objective should use stdout or an seperate outfile, default: False 
            - objargs   (list, optional): List of arguments passed to the objective class at init, default: [] 
            - objkwargs (dict, optional): List of keyword arguments passed to the objective class at init, default: {} 
            - mmkwargs  (dict, optional): List of keyword arguments passed to mm machine at setup, default: {}, in this case the
                defaults in defmmkwargs are used
        """
        # in this dictionary the implemeted objectives has to be given
        objclasses = {'ric_fit': ric_fit.ric_fit, 'force_ric_fit': force_ric_fit.force_ric_fit}
        defmmkwargs = {'lammps': {'screen':False,'logfile':'none'},'pydlpoly':{}}
        assert self.bsetup == False, "FFgen already set up. Not possible to add more objectives"
        assert objname in objclasses.keys(), "Objective %s not available" % objname
        assert mm in mmclasses, 'Requested MM backend %s not available' % mm
        if objname == 'force_ric_fit' and self.mpi_size>1: raise NotImplementedError('No MPI fitting with force_ric_fit')
        dir_prefix = {
                'ric_fit':'rf',
                'force_ric_fit':'frf'}
        # do MPI_stuff
        if self.local_size > 1:
            obj_comm = self.local_comm.Split(self.local_rank, 0)
        else:
            obj_comm = self.local_comm
        obj_size = obj_comm.Get_size()
        obj_rank = obj_comm.Get_rank()
        # pass input
        if coord == None: coord = '%s.mfpx' % sysname
        if fpar == None: fpar = sysname
        if reffile == None: reffile = '%s.hdf5' % sysname
        # handle filenames
        if self.mpi_rank!=0:
            name = '%s_%s-r%i' % (dir_prefix[objname], sysname, self.mpi_rank) 
        else:
            name = '%s_%s' % (dir_prefix[objname], sysname)
        if stdout:
            outfile = sys.stdout
        else:
            outfile=open('%s.out' % name,'w')
        # init molsys
        m = molsys.mol(mpi_comm =obj_comm, out = outfile)
        m.read(os.path.join(self.start_dir,coord))
        # add ff and distribute params
        m.addon('ff', par = self.par)
        if offline:
            assert refsys is not None 
            m.ff.load_params_from_parfile(os.path.join(self.start_dir,fpar), fit=True)
            m.ff.assign_params_offline(refsys)
        else:
            m.ff.read(os.path.join(self.start_dir,fpar), fit=True) 
        # if the current local rank is not responsible for the current obj, add
        # add None to the list of objs
        if self.local_size >1 and self.local_rank != len(self.objs):
            self.add_objective(None, weight)
            return
        # add calculator
        # introduce global rank specific, objective and system dependent name for creating subdirs
        if len(mmkwargs) == 0: mmkwargs = defmmkwargs[mm]
        calc = mmclasses[mm](name, mpi_comm = obj_comm, out = outfile)
        calc.setup(mol=m, local = False, bcond = bcond, **mmkwargs)
        # init objective
        obj = objclasses[objname](calc, os.path.join(self.start_dir,reffile),tag,self.start_dir,out=outfile, mpi_comm=obj_comm,*objargs, **objkwargs)
        # add objective
        self.add_objective(obj, weight)
        return

    def add_objective(self, obj, weight = 1.0):
        """
        Method to add an objective function to the list of objectives

        :Parameters:
            - obj: objective class instance
            - weight(float, optional): weight of the objective in the overall objective function, default: 1.0
        """
        assert self.bsetup == False
        self.objs.append(obj)
        self.weights.append(weight)
        return

    def setup(self):
        """
        Method to setup the final objective function and weights, after 
        calling this method no more objectives could be added. Total sum
        of individual weights has to be zero
        """
        self.bsetup = True
        self.results = np.zeros([len(self.objs)])
        #self.detailed_results = []
        self.weights = np.array(self.weights)
        self.pmin = self.par.variables.ranges[:,0]
        self.pmax = self.par.variables.ranges[:,1]
        self.valrange = self.pmax - self.pmin
        self.initials = copy.deepcopy(self.par.variables.vals)
        self.bounds = self.get_bounds()
        if np.isclose(np.sum(self.weights),1.0) == False:
            self.weights = [1./len(self.objs) for i in range(len(self.objs))]
        # do mpi stuff
        self.nobjperrank = len(self.objs)/self.local_size
        self.objidx = np.array(range(self.local_rank*self.nobjperrank, 
            (self.local_rank+1)*self.nobjperrank))
        self.rresults = np.zeros([self.nobjperrank])
        self.robjs = [self.objs[i] for i in range(len(self.objs)) if i in self.objidx]
#        import pdb;pdb.set_trace()
        assert len(self.results) == len(self.objs) == len(self.weights)
        if self.mpi_rank == 0:
            for o, w in zip(self.objs, self.weights):
                self.status_print('Global', "Obj %s with weight %4.3f applied" % (o, w))
        return


    def set_logger_level(self, level='INFO'):
        """
        Method to set logger level

        :Parameters:
            - level(string, optional): loggerlevel, defaults to INFO
        """
        if level == 'INFO':
            logger.setLevel(logging.INFO)
        elif level == 'WARNING':
            logger.setLevel(logging.WARNING)
        elif level == 'ERROR':
            logger.setLevel(logging.ERROR)
        elif level == 'DEBUG':
            logger.setLevel(logging.DEBUG)
        return

    @timer("calc objectives")
    def __call__(self, params):
        """
        Method to calculate the final objective function, by its 
        individual contributions

        :Parameters:
            - params        (list): list of parameters

        :Returns:
            - msd           (float): fitness 
        """
        assert self.bsetup == True
        self.par.variables(self.convert_from_range(params))
        for i, obj in enumerate(self.robjs):
            msd = obj()
            if self.minimize:
                self.rresults[i] = msd
            else:
                self.rresults[i] = 1./msd
        self.local_comm.Allgatherv(self.rresults,self.results)
        return np.sum(self.results*self.weights)


    def random_guess(self):
        """
        Set up an initial random guess of the parameters, 
        within the parameter ranges

        :Returns:
            - params        (list): list of parameters
        """
        params = []
        for i in range(len(self.pmin)):
            params.append(np.random.uniform(self.pmin[i], self.pmax[i]))
        return self.convert_to_range(params)

    def half_guess(self):
        """
        Set up an initial guess of the parameters, which lys in the half 
        of the individual ranges
        
        :Returns:
            - params        (list): list of parameters
        """
        params = []
        for i in range(len(self.pmin)):
            params.append(self.pmin[i]+0.5*self.valrange[i])
        return self.convert_to_range(params)

    def initial_guess(self):
        """
        Set up an initial guess based on the initial values
        
        :Returns:
            - params        (list): list of parameters
        """
        return self.convert_to_range(self.initials)

    def convert_to_range(self, params):
        """
        Convertes a set of unscaled parameters into a set of parameters with values between 0 and 1
        
        :Parameters:
            - params        (list): list of unscaled parameters

        :Returns:
            - params        (list): list of scaled parameters
        """
        if self.use_range:
            params = (params - self.pmin)/self.valrange
        return params

    def convert_from_range(self,params):
        """
        Convertes a set of scaled parameters to their actual values

        :Parameters:
            - params        (list): list of scaled parameters

        :Returns:
            - params        (list): list of unscaled parameters
        """
        if self.use_range:
            params = self.pmin +(self.valrange*params)
        return params

    def get_bounds(self):
        """
        Method to extract the bounds for the parametersets from
        the ff addon instance of molsys. Here the strings "i" and "h" and "z"
        are translated to weak or hard bounds.

        :Returns:
            - list of tuples, containing the boundary values
        """
        bounds = []
        if self.use_range:
            pmin = self.convert_to_range(self.pmin)
            pmax = self.convert_to_range(self.pmax)
        else:
            pmin = self.pmin
            pmax = self.pmax
        for i,v in enumerate(self.par.variables.values()):
            if v.bounds == ["i","i"]:
                bounds.append([-np.inf, np.inf])
            elif v.bounds == ["h", "i"]:
                bounds.append([pmin[i], np.inf])
            elif v.bounds == ["i", "h"]:
                bounds.append([-np.inf,pmax[i]])
            elif v.bounds == ['h','h']:
                bounds.append([pmin[i],pmax[i]])
            elif v.bounds == ['z','h']:
                bounds.append([0.0,pmax[i]])
            elif v.bounds == ['z','i']:
                bounds.append([0.0,np.inf])
        return zip(*bounds)

    @timer("finish")
    def finish(self, params):
        """
        Method to call after the optimization is finished to write out
        the optimized parameters and further statistics, defined by the
        used objective module.

        :Parameters:
            - params (list): list of optimized parameters
        """
        self.__call__(params)
        self.objs[0].pd.mol.ff.write('opt')
        for obj in self.robjs:
            obj.finish()
        self.finished = True
        return

    def add_validation(self, dirname):
        """
        Method to add an external validations. This external validations have to be defined in an
        own directory. After the fitting has been finished, this directory is copied to rundir. In
        addition the opt.fpar and opt.ric files are copied to the validation directory. The validation
        directory has to contain a file called validate.py.

        :Parameters:
            - dirname (str): Name of the validation directory which is added
        """
        assert type(dirname) == str
        path = os.path.join(self.start_dir, dirname)
        assert os.path.isdir(path) == True, 'Validation directory %s not found' % dirname
        assert os.path.isfile(os.path.join(path,'validate.py')) == True, 'validate.py script not found'
        self.validations.append(dirname)
        return

    @timer('validate')
    def validate(self):
        """
        Method to execute the predefinded external validations.
        """
        assert self.finished == True
        if self.mpi_rank != 0: return
        # run over all validations
        for i,v in enumerate(self.validations):
            path = os.path.join(self.start_dir, v)
            valdir = string.join(self.rundir,v)
            # copy directory to rundir
            shutil.copytree(path, self.rundir)
            # copy opt.fpar and opt.ric to the new directory
            shutil.copy(string.join(self.rundir,'opt.fpar'), valdir)
            shutil.copy(string.join(self.rundir,'opt.ric'), valdir)
            # change into the valdir
            os.chdir(valdir)
            # execute the validation script
            subprocess.call(['python', 'validate.py'], stdout='validation.out', stderr='validation.err')
            # change back to rundir
            os.chdir(self.rundir)
        return
