# -*- coding: utf-8 -*-
'''
                     pydlpoly

(C) 2010-2014 CMC group, Rochus Schmid, Ruhr-Uni Bochum

This is the main frontend to the Python wrapped DL_Poly (V2): PyDLP
It uses the F2PY generated shared object library _pydlpoly, which contains
the DL_Poly code as well as some addtional fortran routines for data
access.
This is a revised version converted into a class
set up is done in init and we include also the conversion from the default
input files. note that we can instantiate only one object because of the
global fortran data structures

revisions:
    - the inner fortran energy routine dlp.calc_energy_force() is only
      called from one single point (calc_energy).
      if we want to get forces as numpy array we can call calc_energy_force, which
      also calls clac_energy. this means also to remove the run() method, which calls
      the dlpoly mainloop and bypasses therefore calc_energy.
      as a consequence all changes to the energy function on python level (QEq, MetaMD etc)
      can be hooked into calc_energy and will ALWAYS be called! (RS, Sept 2011)
    - Instead of the original MPI imp we use now mpi4py as a harness. This means MPI_Init is called
      upon importing MPI from mpi4py. This allows pure python level MPI operations without recurse to the
      Fortran routines in basic_comms.f used by the dl_poly routines.
    - the hdf5 trajectory part is now removed and put into the external class pdlpio
      this class can not just store trajectory data but also restart data and system information which can be used to restart the
      calculation without a tinker xyz file  (October 2012)
    - added MEXT by Tod A. Smith to allow running multiple instances of pydlpoly

NOTE: the molecules features should at some point be seperated into a
class for molecules

'''
import copy
import numpy
import numpy.random as nrand

from mpi4py import MPI

numpy.set_printoptions(threshold=2000)

import re
import string
import time

import os
import os.path
import sys
import shutil
import types
import h5py

import timer

import mext
import ff2pydlpoly

from molsys import mpiobject

# check if we have plumed available
try:
    import _pyplumed
except ImportError:
    plumed_available = False
else:
    plumed_available = True

## import DL_Poly
#import _pydlpoly
## map module names to shortcuts
#dlp       = _pydlpoly.pydlpoly
#dlp_conf  = _pydlpoly.config_module
#dlp_setup = _pydlpoly.setup_module
#dlp_excl  = _pydlpoly.exclude_module
#dlp_vdw   = _pydlpoly.vdw_module
#dlp_time  = _pydlpoly.timing_module
#dlp_prop  = _pydlpoly.property_module
#dlp_mol   = _pydlpoly.molecule_module
#dlp_coul  = _pydlpoly.coulomb_module
#dlp_comm  = _pydlpoly.comm_module
#dlp_nlst  = _pydlpoly.nlist_builders_module
#dlp_efld  = _pydlpoly.external_field_module
#dlp_bnd   = _pydlpoly.bonds_module
#dlp_ang   = _pydlpoly.angles_module
#dlp_dih   = _pydlpoly.dihedral_module
#dlp_inv   = _pydlpoly.inversion_module

#pydlpoly_initalized = False

# further stuff needed
import QEq as qeq
import assign_FF
import unit_cell
import lbfgs
from vectortools import rotate_random

# class for molecule objects
import pdlpmol
# class for IO
import pdlpio

import elements

# dl_poly integer values for different ensembles
#    for all appart "nve" we need a thermostat ("ber" or "hoover", add 1 for "hoover")
enstype = { "nve" : 0, "nvt" : 2, "npt" : 4, "nst" : 6 }

# prsunt is the conversion for pressure from internal units to katm as defined in setup_module.f
prsunt=0.163882576
# boltzmann factor as defined in setup_moldule.f
boltz=8.31451115e-1
# boltzmann factor in kcal/mol/K
boltzkcal = 0.0019872041
# ang2bohr
ang2bohr = 1.0/0.52917721092
# kcal2au
kcal2au = 1.0/627.509
# g/mol2au
g2au = 1.097769e27/6.023e23


class PydlpolyError(Exception):

    def __init__(self, *args, **kwargs):
        Exception.__init__(self,*args,**kwargs)

    #RS: presumably the arguments are storde in self.args (that is what i udnerstand from the docu)
    #    so in to_dict it should be possible to read that out and put it in the dictionary.
    #    NEEDS TESTING

    def to_dict(self):
        rv = {}
        rv["error"]="PydlpolyError"
        rv["message"]="Something went wrong in Pydlpoly"
        rv["code"]=self.args
        return rv



class pydlpoly(mpiobject):

    dlp_instances = 0

    def __init__(self, name, useMExt = False, pseudoRandom = False, mpi_comm = None, out = None):
        ''' startup dlpoly using name as a general runname
            just starting up ... no system loaded, yet.
        '''
        # we use the mext module (and its MExt class) to generate independent instances of pydlpoly
        # Note that the parallelization (as well as possible instances) are identical for each pydlpoly object
        # get _pydlpoly.so module imported as an object self.dlp
        super(pydlpoly,self).__init__(mpi_comm, out)
        if useMExt:
            self._pydlpoly = mext.MExt("_pydlpoly")
            self._pdlp_mext = True
            pydlpoly.dlp_instances += 1
        else:
            assert pydlpoly.dlp_instances == 0, 'No multiple instances without MExt!'
            pydlpoly.dlp_instances = 1
            import _pydlpoly
            self._pydlpoly = _pydlpoly
            self._pdlp_mext = False
        self.dlp       = self._pydlpoly.pydlpoly
        self.dlp_conf  = self._pydlpoly.config_module
        self.dlp_setup = self._pydlpoly.setup_module
        self.dlp_excl  = self._pydlpoly.exclude_module
        self.dlp_vdw   = self._pydlpoly.vdw_module
        self.dlp_time  = self._pydlpoly.timing_module
        self.dlp_prop  = self._pydlpoly.property_module
        self.dlp_mol   = self._pydlpoly.molecule_module
        self.dlp_coul  = self._pydlpoly.coulomb_module
        self.dlp_comm  = self._pydlpoly.comm_module
        self.dlp_nlst  = self._pydlpoly.nlist_builders_module
        self.dlp_efld  = self._pydlpoly.external_field_module
        self.dlp_bnd   = self._pydlpoly.bonds_module
        self.dlp_ang   = self._pydlpoly.angles_module
        self.dlp_dih   = self._pydlpoly.dihedral_module
        self.dlp_inv   = self._pydlpoly.inversion_module
        ######## Attach Error callback function ######################
        self._pydlpoly.pyerror = self.error
        ######## Initialize random number generators #################
        if pseudoRandom:
            if pseudoRandom is True:
                self.pseudo_random()
            else:
                self.pseudo_random(seed=pseudoRandom)
        else:
            self.true_random()
        ######### Now it is safe to start up ##########################
        self.local_comm = self.mpi_comm
        self.dlp_comm.pydlpoly_comm = self.local_comm.py2f()
        self.instances = None
        # set some basics defaults here
        self.eprec = 4
        self.eprec_incr = 4
        # Now we are set to init everything
        self.dlp.init_dlpoly()
        self.idnode = self.dlp.idnode
        self.nodes  = self.dlp.mxnode
        # module level global variables
        self.name = name
        self.natoms = 0
        self.loc_natoms = 0
        #
        self.keep_raw_files = False
        #ccs
        # These two are only needed for a dlp lbfgs qeq calculation...
        self.counter = 0
        self.H = numpy.eye(self.natoms)
        #
        self.qeq = None
        #
        self.QMMM = False
        #
        # extra_term can be used to add an additional energy term to the systems energy expression
        self.extra_term = None
        #
        # extra_systems can be a list of dynamic systems which must provide at least
        # the following methods: calc_energy(), vv_fh(), vv_sh(), get_ekin(), get_epot(), get_name(), get_dof()
        # currently extra systems are only considered during MD and are propagated via vv_fh/sh
        # with a vleocity verlet in parallel to the atoms
        self.extra_systems = []
        #
        self.pdlpio = None
        #
        self.plumed = False
        self.plumed_init = False
        #ccs
        self.set_atoms_moved()
        #
        self.enforce_nlist_rebuild = False
        #
        self.nhosts=0
        #
        self.use_Jij=False
        self.valid_Jij=False
        #
        self.trackimg=False     
        #
	# for constraints
        self.constraint = False
        self.use_efield = False
        #
        self.use_split_vdw = False
        #
        # defaults for the current md simaulation stage (just to keep things in python)
        self.md_counter = 0
        self.md_ens     = "nve"
        self.md_thermo  = None
        self.md_T       = None
        self.md_p       = None
        self.md_relax   = None
        self.md_stage   = None
        self.md_traj    = None
        self.md_rnstep  = 10
        self.md_tnstep  = 100
        self.md_grlmb   = None
        # default settings for CONTROL
        self.control = {\
        "timestep" : 0.001 ,\
        "temp"     : 300.0 ,\
        "cut"      : 12.0  ,\
        "delr"     : 1.0
        }
        # runmode
        self.mode=None
        # virtual atoms
        self.virtual_atoms = None
        #
        self.print_header()
        #
        # init a timer class
        self.timer = timer.timer("pydlpoly")
        #
        if True_Random==True:
            self.pprint("Using true random in python and fortran")
        elif True_Random==False:
            self.pprint("Using pseudo random in python and fortran -> for testing")
        else:
            raise ValueError, "Undefined random state"
        return

    @property
    def is_master(self):
        return self.idnode == 0
        
    def error(self):
        self.pprint("PYDLPOLY ERROR STATE")
        self.pprint("Error code %d" % self.dlp_comm.error_code)
        raise PydlpolyError(self.dlp_comm.error_code)
        return

    def true_random(self):
        global True_Random
        # default behavior of python is true random anyway
        self.dlp.init_dlpoly_random()
        True_Random = True
        return

    def pseudo_random(self,seed=42):
        global True_Random
        # keep pseudo random of fortran but seed numpy random generator
        nrand.seed(seed)
        True_Random = False
        return

    def print_header(self):
        # printout header
        self.mark_line = 80*"*"+"\n"
        self.empty_line = 80*" "+"\n"
        header =  """
            PYDLPOLY

            pydlpoly is a Python wrapped and extended fork of DL_Poly
            (Version 2 now referred to as DL_Poly Classic)
            developed by the CMC group

            Concept and lead development: Rochus Schmid

            Contributions by:
            Christian Spickermann, Sareeya Bureekaew, Mohammad Alaghemandi,
            Johannes Dürholt, Dennis Pache, Julian Keupp, Niklas Siemer
            Copyright Ruhr-Universität Bochum 2010-2015

            Original Copyright for the DL_Poly basis:
            ##################################################################
            dl_poly classic is an stfc/ccp5 program package for the
            dynamical simulation of molecular systems.

            dl_poly is the copyright of the stfc daresbury laboratory,
            daresbury, warrington wa4 4ad.

            neither the stfc, daresbury laboratory, ccp5 nor the authors
            of this package claim that it is free from errors and do not
            accept liability for any loss or damage that may arise from
            its use. it is the users responsibility to verify that the
            package dl_poly is fit for the purpose the user intends for
            it.

            users of this package are recommended to consult the dl_poly
            user manual for the full description of its use and purpose.

            authors: w.smith and t.r.forester 1995
            copyright daresbury laboratory 1995
            ####################################################################
        """
        self.pprint(self.mark_line)
        self.pprint(self.empty_line)
        self.pprint(header)
        self.pprint(self.mark_line)
        self.pprint(self.empty_line)
        self.pprint("\n\nWelcome on board!!!!")
        self.pprint("running on %d nodes\n" % self.dlp.mxnode)
        return

    def get_comm(self):
        return self.local_comm

    ########################################################################
    # various setup actions
    #    raw : use exisiting standard files in the current directory
    def setup_raw(self):
        self.pprint(self.mark_line)
        self.pprint("Running in local directory using standard filenames\n")
        self.pprint(self.mark_line)
        self._startup_dlpoly()
        return

    def setup(self, key=None, control=None, xyz=None, path=None, local=False, \
        molecules=None, bcond=None, QMMM=None, empty_box=None,\
        restart=False, pdlp=None, read_stage="default", start_stage="default", velocities=False,\
        imgidx = False, pdlp_restart=None, do_first_energy=True, smallrings=False, keep=False, vdw_srcut=None,\
        split_vdw=False, warning_FF=False, mol = None, par = None, ff = None):
        """
        Setup method to read input files of the system (or restart) and to initalize pydlpoly
        
        If mol or par or ff is provided then the new setup via molsys is invoked

        :Parameters:
            - key             (str) : filename of key file (if `None` <runname>.key is used)
            - control         (str) : filename of control file (if `None` default is used)
            - xyz             (str) : filename of xyz file (if `None` <runname>.xyz is used)
            - path            (str) :
            - local           (bool): default: False - generate a subdir from runname and run job there (True - run local)
            - molecules       (list): add guest molecules
            - bcond           (int) : enforce boundary conditions (1: cubic, 2: orthorombic, 3: triclinic, 6: xy-periodic)
            - QMMM
            - empty_box
            - restart
            - pdlp
            - read_stage
            - start_stage
            - velocities
            - imgidx
            - pdlp_restart
            - do_first_energy (bool): default True  - compute a first energy after reading the system in
            - smallrings      (bool): default False - skip detecting small rings (4, 5 rings) for e.g. angle5, bond5 types
            - keep            (bool): Keep raw DL_Poly files (CONFIG, CONTROL, FIELD, OUTPUT), default=False
            - vdw_srcut       (float): default is None - if a float is specified this is the cutoff energy in kcal/mol above which vdw rep is replaced by an inverted parabola
            - split_vdw       (bool): splitting vdw into dispersive and repulsive part (for seperate scling in GCMD) default: False
            - warning_FF      (bool): ???? who did that ???
            - mol             (molsys instance): default None
            - par             (string): name of the parameter files .ric/.par
            - ff              (string): name of the force field to be assigned (if ff="file" we use default name for par)
        """
        # check correct input for new setup via molsys
        legacy_mode = True         # this is true if assign_FF is used and False if molsys is used
        if mol is not None: 
            assert par ==  ff == None
            legacy_mode = False
        if par is not None: 
            assert mol ==  ff == None
            legacy_mode = False
        if ff  is not None: 
            assert mol == par == None
            if ff == "file":
                par = self.name
            legacy_mode = False
        # setup proper filenames
        self.local = local
        self.start_dir = os.getcwd()
        if xyz:
            self.xyz_file = xyz
        else:
            if legacy_mode:
                self.xyz_file = self.name + ".xyz"
            else:
                self.xyz_file = self.name + ".mfpx"
        if key:
            self.key_file = key
        else:
            self.key_file = self.name + ".key"
        # use control only if requested otherwise provide defaults from directory self.control
        self.control_file=control
        if path:
            self.path = path
        else:
            self.path = self.start_dir
        if pdlp:
            self.pdlp_file = pdlp
        else:
            self.pdlp_file = self.name + ".pdlp"
        if pdlp_restart:
            self.pdlp_rest_file = pdlp_restart
        else:
            self.pdlp_rest_file = self.pdlp_file
            if molecules and restart:
                raise IOError, "You need to specify a new pdlp file if you add molecules!"
        #
        if not local:
            self.rundir = self.path + "/" + self.name
            if self.local_comm.Get_rank() == 0:
                # check the presence of a rundir
                if os.path.isdir(self.rundir):
                    # directory exists - increment
                    i = 1
                    temprundir = self.rundir + ("_%d" % i)
                    while os.path.isdir(temprundir):
                        i += 1
                        temprundir = self.rundir + ("_%d" % i)
                    self.rundir = temprundir
                # now we know that this directory does not exist and can make it
                os.mkdir(self.rundir)
            # at this point the name of the actual rundir is known only on the GLOBAL master
            # thus we broadcast it to the nodes before we cd
            self.rundir = self.local_comm.bcast(self.rundir)
            # we need to sync here in order to let the master, making the directory, catch up
            self.local_comm.Barrier()
            os.chdir(self.rundir)
        else:
            # we want to run in the local directory - so basically we do not need to do anything
            # NOTE we need to check if multiple instances are used. in this case it is not safe to use local
            if pydlpoly.dlp_instances > 1:
                raise IOError, "This is pydlpoly instance %d, using local=True in setup is not safe!" % dlp_instances
            self.rundir=self.start_dir
        ## write out mode and file information
        self.pprint(self.mark_line)
        if not restart:
            if empty_box:
                if not molecules:
                    raise IOError, "Setting up from empty box without adding molecules leads to empty system! Aborting!"

                self.pprint("Setting up from empty box and added molecules\n")
            else:
                if legacy_mode:
                    self.pprint("Standard setup in legacsy mode from xyz/key file\n")
                else:
                    self.pprint("New Setup using molsys object")
        else:
            self.pprint("Restart using pdlp file\n")
        self.pprint(self.mark_line)
        if restart:
            self.pprint("pdlp file    : %s" % self.pdlp_rest_file)
            self.pprint("       restart data from stage %s" % read_stage)
            if velocities:
                self.pprint("       reading also velocities in")
            if imgidx:
                self.pprint("       reading also image info in")
        else:
            if empty_box:
                self.pprint("box size     : %s" % str(empty_box))
            else:
                self.pprint("xyz file     : %s" % self.xyz_file)
        self.pprint("key file     : %s" % self.key_file)
        if self.control_file:
            self.pprint("control file : %s" % self.control_file)
        else:
            self.pprint("control file : generated from defaults")
            # Select shift damping as default if no other option is present
            elecopts = ['shift', 'spme', 'ewald', 'hke']
            if all(elec not in self.control for elec in elecopts):
                self.control['shift'] = 'damping'
        self.pprint("running in   : %s" % self.rundir)
        if keep:
            self.keep_raw_files = True
            self.pprint("raw DL_Poly files will not be removed")
        self.pprint("\n")
        ## start up and write CONFIG, FIELD and CONTROL
        #
        #  We have two modes: either restart=False and we read
        #       from xyz as usual. in this case also molecules can be added.
        #       if restart=True we read from a pdlp file and we can NOT add
        #       further molecules. In this case we open the pdlp file already here
        #
        if mol is not None:
            self.mol = ff2pydlpoly.wrapper(mol)
        elif par is not None or ff is not None:
            import molsys
            print self.xyz_file
            mol = molsys.mol.fromFile(os.path.join(self.start_dir,self.xyz_file))
            mol.addon('ff')
            if par:
                mol.ff.read(os.path.join(self.start_dir,par))
            else:
                mol.ff.assign_params(ff)
            self.mol = ff2pydlpoly.wrapper(mol)
        else:
            self.mol = assign_FF.mol(self.is_master)
        if not restart:
            if empty_box:
                self.mol.empty_box(empty_box)
            else:
            # read in from tinker file
                self.mol.read_tinker_xyz(self.start_dir+"/"+self.xyz_file)
        else:
            # read from pdlp file ... so open it here already ... we read
            # on all nodes in this case
            pdlpf  = pdlpio.pdlpio(self.start_dir+"/"+self.pdlp_rest_file, mode="r")
            self.mol.get_system_from_pdlp(pdlpf)
            # now we need to get the structure from the restart section
            self.mol.read_pdlp_xyz(pdlpf, read_stage, velocities, imgidx)
            pdlpf.close()
        ## if  is true we expect a list of tuples with (name, file, N)
        added_molnames = []
        if molecules:
            if velocities:
                self.pprint("no velocities if molecules are added!")
                velocitites = False
            for m in molecules:
                mname = m[0]
                mfile = m[1]
                mN    = m[2]
                if len(m) > 3:
                    offset = numpy.array(m[3])
                    scale  = numpy.array(m[4])
                else:
                    offset = numpy.zeros([3], dtype="float64")
                    scale  = numpy.ones([3], dtype="float64")
                self.mol.add_mol_tinker_xyz(self.start_dir+"/"+mfile, mN, mname, offset, scale)
                added_molnames.append(mname)
        # find connectivity and potential terms
        self.mol.load_FF(self.start_dir+"/"+self.key_file)
        self.mol.find_internals(do_smallring_check=smallrings, do_lin_check = False)
        #self.mol.load_FF(self.start_dir+"/"+self.key_file)
        self.mol.assign_FF(warning_FF)
        # assing molecule names (if there are any)
        if not restart:
            self.mol.assign_molnames()
        # if there are QMatoms remove them from the molecule they belong to and make a new one
        if QMMM:
            QMatoms, self.QMMM_interface, self.QMMM_links = QMMM
            self.QMMM = True
            QMatoms = (numpy.array(QMatoms)-1).tolist()
            self.QMmol = self.mol.make_extra_mol(QMatoms, "qm")
            self.mol.exclude_mol("qm")
            # set up the QMMM interface
            self.QMMM_interface.init(self, self.QMmol, self.QMMM_links)
        # if there are rigid atoms exclude them
        self.mol.exclude_rigid()
        # if there are frozen molecules freeze and exclude them
        self.mol.exclude_frozen()
        # print report
        self.mol.report_molecules()
        # do the file writing only on the master node
        if self.is_master:
            self.mol.write_FIELD()
            self.mol.write_CONFIG(bcond=bcond, vel=velocities)
            if self.control_file:
                shutil.copy(self.start_dir+"/"+self.control_file, "CONTROL")
            else:
                self.write_default_control()
        # to make sure that files are written properly from the master and all nodes can startup we wait here for all to catch up
        self.local_comm.Barrier()
        # DEBIUG DEBUG EXPERIMENTAL
        if vdw_srcut:
            # this is a hack .. at this point dlp.engunit is not yet defined .. we need to imply kcal/mol
            self.dlp_vdw.srdamp_ener = vdw_srcut*418.4
            # dlp_vdw.srdamp_ener = vdw_srcut*dlp.engunit
            self.dlp_vdw.lsrdamp = True
            self.pprint("Using a short range damped repulsion for vdw truncating to a inverted parabola above %10.5f kcal/mol" % vdw_srcut)
        if split_vdw:
            self.dlp_vdw.lspltvdw = True
            self.use_split_vdw = True
            self.pprint("Using a split vdw potential. Dispersion and repulsion are computed seperately. This is useful only for GCMD!")
        # DEBUG DEBUG
        #########################################################################################
        #
        #    ok .. done lets start pydlpoly on all nodes
        self._startup_dlpoly(allocate=True)
        #########################################################################################
        #
        #   Molecules startup (using pdlpmol)
        # NOTE: self.added_mols is now a list of molecule objects (and not a list of integers)
        #
        self.molecules = pdlpmol.molecules(self)
        # if there have been molecules added to the system (randomly) we switch them off
        if len(added_molnames)>0:
            self.pprint("\nAdded molecules will be set to lambda=zero (switched off)")
            self.added_mols = []
            for mname in added_molnames:
                mollist = self.molecules.get_list_from_name(mname)
                self.added_mols += mollist
                self.pprint("%10s: %d" % (mname, len(mollist)))
                if self.use_split_vdw:
                    for m in mollist: m.set_lambda(0.0, 0.0, 0.0)
                else:
                    for m in mollist: m.set_lambda(0.0, 0.0)
        #####     initalize virtual atoms if there are any (needs to be done after statup of fortan engine
        # RS NOTE  should at some point be done in Fortran (or Cython??)
        if len(self.mol.virtual_atoms) > 0:
            self.pprint("Initalizing %d virtual atom sites" % len(self.mol.virtual_atoms))
            self.virtual_atoms = self.mol.virtual_atoms
            self.virt_atom_defs = self.mol.virt_atom_defs
            self.frozen_virtuals = len(self.virtual_atoms)*[False]
            self.init_virt_atoms()
        #####      QMMM startup of molecules (needs to be done after init of molecules
        if self.QMMM:
            # this is a tricky one here: now we have molecules initalized
            # this means we can set lambdas ... here we switch off the internal coul and vdw interactions
            self.dlp_mol.mol_lambint_vdw[self.QMmol]  = 0.0
            self.dlp_mol.mol_lambint_coul[self.QMmol] = 0.0
            self.QMMM_interface.setup()
        # at this point we write the system data from the mol instance to the restart file (only on master)
        self.pdlpio = False
        self.pprint("Writing Trajectory and Restart data to pdlp file (hdf5 format)")
        self.pprint("File %s is written in rundir %s" % (self.pdlp_file, self.rundir))
        if self.is_master:
            self.pdlpio = pdlpio.pdlpio(self.pdlp_file)
            self.mol.write_system_to_pdlp(self.pdlpio)
            # set up datasets and groups in the default stage
            self.pdlpio.initialise(pd=self, initial_stage=start_stage)
            # FOR DEBUG
            # self.pdlpio.verbose = 1
            # write intitial geometry to the restart section
            self.pdlpio()
        # compute first energy
        # MExt ... go back to startdir .  this is needed in order to have a predictable directory when multiple instances of pydlpoly run
        os.chdir(self.start_dir)
        if imgidx: self.track_images(imgidx= self.mol.imgidx)
        if do_first_energy:
            self.calc_energy()
            self.pprint(self.mark_line)
            self.pprint("Energy of initial configuration computed")
            self.report_energies(full=True)
        return

    def write_default_control(self):
        """ writes a default control file to the current directory.
            you can change settings by modifying the directory.
            each key is written if it is not None.
            So to switch off a setting set it to None.
            This is for convenience and some things can not be changed ... use your own CONTOL if you want to override defaults
            NOTE: The temperature setting is actually meaningless because the temperature of an MD simulation is set with MD_init!

            some keywords want to know the cutoff to be set already (e.g. spme) .. so we enforce writing it first
        """
        if self.is_master:
            f = open("CONTROL", "w")
            f.write("SYSTEM: %s\n" % (self.name.upper()))
            f.write("\n")
            keys = self.control.keys()
            # remove cut and delr and add it to the front ... they are defaults and always present
            keys.remove("cut")
            keys.remove("delr")
            keys = ["cut", "delr"]+keys
            for k in keys:
                if self.control[k] :
                    f.write("%20s  %s\n" % (k, str(self.control[k])))
            f.write("\n")
            # these checks are check if we have a MOF-FF type FF
            # check if we have gaussian charges in the key file
            if self.mol.gaussians_used:
                f.write("%20s\n" % "gaussian")
                if self.mol.spbasis_used:
                    f.write("%20s\n" % "sp-basis")
#                if self.mol.FF.settings["chg-13-scale"] == 1.0 and self.mol.FF.settings["chg-12-scale"] == 1.0 :
#                    f.write("%20s\n" % "all cterms")
            # now write the fixed part
            f.write("%20s\n" % "all cterms")
            f.write(\
                """
                sw_vdw  0.9
                sw_chrg 0.9
                derivative vdw spline

                finish

                """)
            f.close()
        return


    def _startup_dlpoly(self, allocate=True):
        """ This should not be used directly -> use via setup strategies """
        # now we are in the right directory and CONTROL, FIELD and CONFIG are present
        if allocate: self.dlp.allocate_dlpoly()
        self.dlp.setup_dlpoly()
        self.natoms = self.dlp_setup.mxatms
        self.pprint(self.mark_line)
        self.pprint("pydlpoly is starting up!\nsystem with %d atoms loaded" % self.natoms)
        # enforce all optimizer switches to be zero
        self.dlp.loptim  = False
        self.dlp.lzero   = False
        self.dlp.lminim  = False
        self.dlp_time.timer_init()
        self.pprint("Timers initialized")
        return

    def add_stage(self, stage, traj=None, rnstep=10, tnstep=100):
        if self.is_master:
            if self.pdlpio.has_stage(stage):
                raise IOError, "The PDLP restart file already has a stage of this name!"
            self.pdlpio.add_stage(stage, rnstep)
            if traj:
                self.pdlpio.add_traj(tnstep, traj)
        return

    def reFIELD(self, filename):
        """ This method reloads a system from a given tinker (!!) xyz file
            no molecules, restarts or other filetypes are suppported, only a sinlge tinker xyz file
            in addtion, the xyz file must be an isomer to the inital one (same number of atoms, bonds etc)
            BUT: the order of the atoms can change and also the specific force field terms
            in other words: the xyz and force field info is overwritten """
        self.pprint("\nDoing a reFIELD: rereading system from %s" % filename)
        self.pprint("   WARNING! This implies that the xyz file (Tinker format) must be")
        self.pprint("            an isomer of the already read system (same number of atoms/bonds etc.)")
        self.pprint("            No reallocation of arrays is done but just a new FIELD file is reparsed")
        self.mol.read_tinker_xyz(self.start_dir+"/"+filename)
        self.mol.find_internals()
        self.mol.assign_FF()
        if self.is_master:
            self.mol.write_FIELD()
        self.dlp.reload_field()
        xyz = numpy.array(self.mol.xyz)
        self.set_xyz(numpy.array(self.mol.xyz))
        self.set_cell(self.mol.cell)
        #RS need to quench all  here for MD !!!!
        self.set_atoms_moved()
        return


    #########################################################################
    ## "default" run ... just start the standard dlpoly mainloop
    ##  this function will eventually go!!!!

    def orig_run(self):
        self.pprint("Running dlpoly mainloop ... this is depreciated feature!")
        self.dlp.init_dyn(skip_REVOLD=False)
        self.dlp.dynamics()
        self.pprint("done")
        return

    def end(self):
        self.pprint("your ride with pydlpoly ends here! Thanks for flying with us.")
        # self.report_timing()
        if self.qeq: self.qeq.timer.report("#q#")
        self.timer.report("#p#")
        if not self.keep_raw_files:
            self.pprint("Removing raw DL_Poly input files in rundir (to avoid use 'keep' in setup)")
            if self.is_master:
                os.remove(self.rundir+"/OUTPUT")
                os.remove(self.rundir+"/CONFIG")
                os.remove(self.rundir+"/FIELD")
                os.remove(self.rundir+"/CONTROL")
        #dlp.finish_dlpoly()
        del(self._pydlpoly)
        return

    def __del__(self):
        self.end()
        return

    ######################################################################
    #                core calculation routines
    #  all methods should rely on these to have a central point for
    #  changing things
    ######################################################################

    def calc_energy(self, force=True):
        """ compute the current energy, by default also the force is computed
            Note that in dlpoly the force is always computed (even with force=False)
        """
        if self.atoms_moved:
            self.timer.switch_on("energy")
            if self.qeq:
                self.timer.switch_to("qeq")
                self.qeq.preopt_q()
                self.timer.switch_to("energy")
            if self.virtual_atoms:
                self.timer.switch_to("virt atoms set pos")
                self.virt_atom_set_pos()
                self.timer.switch_to("energy")
            if self.enforce_nlist_rebuild:
                # if some operation (move atoms/molecules around) could break the neigborlist enforce its rebuild
                self.dlp_nlst.nlst_notest  = True
                self.dlp_nlst.nlst_rebuild = True
            ### nullify stress tensor
            self.dlp_conf.stress[:]=0.0
            ### MAIN ROUTINE TO COMPUTE ENERGY AND FORCE IN DL_POLY ###
            self.dlp.calc_energy_force()
            ### ################################################### ###
            ### MAIN ROUTINE TO COMPUTE ENERGY AND FORCE IN DL_POLY ###
            if self.QMMM:
                self.timer.switch_to("qmmm")
                self.QMMM_eng = self.QMMM_interface(force)
                self.timer.switch_to("energy")
                # for a QMMM run we report the energies every step
                self.report_energies()
            if self.virtual_atoms:
                self.timer.switch_to("virt atoms dist force")
                self.virt_atom_distribute_force()
                self.timer.switch_to("energy")
            self.atoms_moved = False
            if self.enforce_nlist_rebuild:
                self.dlp_nlst.nlst_notest  = False
                self.dlp_nlst.nlst_rebuild = False
                self.enforce_nlist_rebuild = False
#            if self.QMMM:
#                self.timer.switch_to("qmmm")
#                self.QMMM_eng = self.QMMM_interface(force)
#                self.timer.switch_to("energy")
#                # for a QMMM run we report the energies every step
#                self.report_energies()
            if self.extra_term:
                self.timer.switch_to("extra term")
                self.extra_eng = self.extra_term()
                self.timer.switch_to("energy")
            if self.plumed:
                self.timer.switch_to("metamd plumed")
                self.meta_eng = self.plumed_energy()
                self.timer.switch_to("energy")
            if self.constraint:
                self.multiply_force(self.constraint_mask)
            for s in self.extra_systems:
                self.timer.switch_to(s.get_name())
                s.calc_energy()
                self.timer.switch_to("energy")
            self.timer.switch_off()
        # return configuration energy in current units
        energy = self.dlp.engcfg/self.dlp.engunit
        if self.QMMM:       energy += self.QMMM_eng
        if self.extra_term: energy += self.extra_eng
        if self.plumed:     energy += self.meta_eng
        return energy

    def rebuild_nlist(self):
        self.enforce_nlist_rebuild = True
        return

    def calc_energy_force(self):
        """
        Computes energy with calc_energy and returns also current forces

        :Returns:

             - energy (energy in kcal/mol)
             - fxyz  [natoms, 3]
        """
        energy = self.calc_energy(force=True)
        fxyz = self.get_force()
        return energy, fxyz

    def get_force(self):
        # for convenience we remap the linear coordinates in a [N,3] shape
        # this could be slow for larger systems ...
        fxyz = numpy.empty([self.natoms, 3],"d")
        fxyz[:,0] = self.dlp_conf.fxx
        fxyz[:,1] = self.dlp_conf.fyy
        fxyz[:,2] = self.dlp_conf.fzz
        fxyz *= 1.0/self.dlp.engunit
        return fxyz

    def get_subset_force(self, aind):
        # for convenience we remap the linear coordinates in a [N,3] shape
        # this could be slow for larger systems ...
        fxyz = numpy.empty([len(aind), 3],"d")
        fxyz[:,0] = self.dlp_conf.fxx[aind]
        fxyz[:,1] = self.dlp_conf.fyy[aind]
        fxyz[:,2] = self.dlp_conf.fzz[aind]
        fxyz *= 1.0/self.dlp.engunit
        return fxyz

    def add_force(self, fxyz):
        """
        force is converted back to internal units and added to the fortran arrays

        :Paramters:
            - fxyz : force array to be added
        """
        self.dlp_conf.fxx[::] += fxyz[:,0]*self.dlp.engunit
        self.dlp_conf.fyy[::] += fxyz[:,1]*self.dlp.engunit
        self.dlp_conf.fzz[::] += fxyz[:,2]*self.dlp.engunit
        return

    def add_subset_force(self, fxyz, aind):
        """
        force is converted back to internal units and added to the fortran arrays

        :Paramters:
            - fxyz : force array to be added
        """
        self.dlp_conf.fxx[aind] += fxyz[:,0]*self.dlp.engunit
        self.dlp_conf.fyy[aind] += fxyz[:,1]*self.dlp.engunit
        self.dlp_conf.fzz[aind] += fxyz[:,2]*self.dlp.engunit
        return

    def multiply_force(self, fxyz):
        self.dlp_conf.fxx[::] *= fxyz[:,0]
        self.dlp_conf.fyy[::] *= fxyz[:,1]
        self.dlp_conf.fzz[::] *= fxyz[:,2]
        return

    def set_force(self, fxyz):
        """
        force is converted back to internal units and written into the fortran arrays

        :Parameters:
            - fxyz : force to be set
        """
        self.dlp_conf.fxx[::] = fxyz[:,0]*self.dlp.engunit
        self.dlp_conf.fyy[::] = fxyz[:,1]*self.dlp.engunit
        self.dlp_conf.fzz[::] = fxyz[:,2]*self.dlp.engunit
        return

    def set_subset_force(self, fxyz, aind):
        """
        force is converted back to internal units and written into the fortran arrays

        :Parameters:
            - fxyz : force to be set
            - aind : array of atom indices to operate on
        """
        self.dlp_conf.fxx[aind] = fxyz[:,0]*self.dlp.engunit
        self.dlp_conf.fyy[aind] = fxyz[:,1]*self.dlp.engunit
        self.dlp_conf.fzz[aind] = fxyz[:,2]*self.dlp.engunit
        return



    def calc_numforce(self,delta=0.0001):
        self.mode=None
        xyz = self.get_xyz()
        energy, fxyz = self.calc_energy_force()
        num_fxyz = numpy.zeros([self.natoms,3],"float64")
        atypes = self.get_atomtypes()
        for a in xrange(self.natoms):
            for i in xrange(3):
                xyz[a,i] += delta
                self.set_xyz(xyz)
                ep = self.calc_energy(force=False)
                ep_contrib = self.get_energy_contribs()
                xyz[a,i] -= 2*delta
                self.set_xyz(xyz)
                em = self.calc_energy(force=False)
                em_contrib = self.get_energy_contribs()
                xyz[a,i] += delta
                num_fxyz[a,i] = -(ep-em)/(2.0*delta)
                # self.pprint("ep em delta_e:  %20.15f %20.15f %20.15f " % (ep, em, ep-em))
                self.pprint(numpy.array2string((em_contrib-ep_contrib)/(2.0*delta),precision=2,suppress_small=True))
                self.pprint("atom %8d (%20s) %3d: anal: %12.6f num: %12.6f diff: %12.6f " % (a,atypes[a],i,fxyz[a,i],num_fxyz[a,i],( fxyz[a,i]-num_fxyz[a,i])))
                os.sys.stdout.flush()
        return fxyz, num_fxyz

    def calc_coul_energy(self, get_core_pot=False, constrain=False, con_mol=False):
        # this will need some improvement in the future ...
        # run just the coulomb part of the energy ... directly via the fortran wrapper to molecular_dynamics
        # fortran call sequence:
        #     molecular_dynamics -> force_manager
        # in force_manager the qeq_iter flag makes it to bypass everything unrelated to Coulomb
        if self.use_Jij:
            # we use Jij storage ... do we have a valid Jij already
            if self.valid_Jij:
                # ok ... just apply
                self.dlp_coul.apply_Jij(self.idnode,self.nodes,self.natoms,self.dlp.loglnk, False, False)
                qpot_fast_Jij = numpy.array(self.dlp_conf.pot[:])/self.dlp.engunit
                #qpot_fast_Jij = self.gdsum(qpot_fast_Jij)
                Ecoul = 0.5*numpy.sum(qpot_fast_Jij*self.dlp_conf.chge)
                return Ecoul,qpot_fast_Jij
            else:
                # no ... prepare for use and do the regular energy call to set up Jij
                self.dlp_coul.zero_Jij()
                self.dlp_coul.store_Jij = True
        self.dlp_setup.qeq_coulonly = True
        #dlp_setup.qeq_iter = True
        self.dlp_setup.calc_hess = False
        if get_core_pot: self.dlp_conf.get_core_pot = True
        self.dlp.calc_energy_force()
        qpot = numpy.array(self.dlp_conf.pot[:])/self.dlp.engunit
        Ecoul = self.dlp.engcpe/self.dlp.engunit
        self.dlp_setup.qeq_coulonly = False
        #dlp_setup.qeq_iter = False
        if self.use_Jij:
            self.dlp_coul.store_Jij=False
            self.valid_Jij=True
        if get_core_pot:
            self.dlp_conf.get_core_pot = False
            core_qpot = numpy.array(self.dlp_conf.core_pot[:])/self.dlp.engunit
            return Ecoul, qpot, core_qpot
        else:
            return Ecoul, qpot

    def calc_Jmat(self, get_core_pot=False):
        self.dlp_setup.qeq_coulonly = True
        self.dlp_setup.calc_hess = True
        if get_core_pot: self.dlp_conf.get_core_pot = True
        self.dlp.calc_energy_force()
        Jmat = numpy.array(self.dlp_conf.q_hess[:,:])/self.dlp.engunit
        qpot = numpy.array(self.dlp_conf.pot[:])/self.dlp.engunit
        Ecoul = self.dlp.engcpe/self.dlp.engunit
        # in the parallel case this needs to be broadcasted??
        self.dlp_setup.qeq_coulonly = False
        self.dlp_setup.calc_hess = False
        if get_core_pot:
            self.dlp_conf.get_core_pot = False
            core_qpot = numpy.array(self.dlp_conf.core_pot[:])/self.dlp.engunit
            return Ecoul, qpot, Jmat, core_qpot
        else:
            return Ecoul, qpot, Jmat

    def get_val_energy(self):
        return self.dlp_conf.enc_val/self.dlp.engunit

    def get_core_energy(self):
        return self.dlp_conf.enc_core/self.dlp.engunit

    def get_coreval_energy(self):
        return self.dlp_conf.enc_coreval/self.dlp.engunit

    def set_extra_term(self, efunc):
        """
        add an additional energy term (called AFTER standard pydlpoly terms)

        :Parameters:
            - efunc : callable object (without paramters) computes an extra energy and forces

        :Note:
            There can be only one extra term
        """
        self.pprint("Extra energy term was added. Recalculating energy")
        self.extra_term = efunc
        self.set_atoms_moved()
        self.calc_energy()
        self.report_energies(full=True)
        return

    def add_extra_system(self, system):
        """
        add an extra system (class instance with a defined API, see docu for extra systems)
        There can be multiple extra systems. Extra systems have there own degrees of freedom which are
        propagated (which is the difference to an extra_term)

        :Parameters:
            - system : instance of a class with the extra system API

        :Note:
            Energy and forces of the extra systen are calculated after the regular pydlpoly energies are computed
        """
        self.pprint("Extra dynamic system %s added" % system.get_name())
        self.extra_systems.append(system)
        system.start_up()
        self.set_atoms_moved()
        return

    def set_atoms_frozen(self, atlist):
        for a in atlist:
            ai = a-1
            self.dlp_conf.lstfrz[ai] = 1
            self.pprint("   == atom %5d frozen (%s, %s)" % (a, self.mol.elems[ai], self.mol.types[ai]))
            self.set_constraint(ai,0)
            self.set_constraint(ai,1)
            self.set_constraint(ai,2)
        return

    def set_atoms_unfrozen(self, atlist):
        for a in atlist:
            ai = a-1
            self.dlp_conf.lstfrz[ai] = 0
            self.pprint("   == atom %5d unfrozen (%s, %s)" % (a, self.mol.elems[ai], self.mol.types[ai]))
            self.release_constraint(ai,0)
            self.release_constraint(ai,1)
            self.release_constraint(ai,2)
        return

    def set_virtuals_frozen(self, atlist):
        if self.virtual_atoms != None:
            for a in atlist:
                ai = a-1
                if ai in self.virtual_atoms:
                    self.frozen_virtuals[self.virtual_atoms.index(ai)] = True
                    self.pprint("  == virtual %5d frozen" % a)
                else:
                    self.pprint("atom %5d is not a virtual site" % (a))
        else:
            self.pprint("There are no virtual sites defined")
        return

    def set_virtuals_unfrozen(self, atlist):
        if self.virtual_atoms != None:
            for a in atlist:
                ai = a-1
                if ai in self.virtual_atoms:
                    self.frozen_virtuals[self.virtual_atoms.index(ai)] = False
                    self.pprint("  == virtual %5d unfrozen" % a)
                else:
                    self.pprint("atom %5d is not a virtual site" % (a))
        else:
            self.pprint("There are no virtual sites defined")
        return

    def set_efield(self, field_vect, const_D=False, use_ref = True, ref = None, reduced = False):
        """
        Switch on finite field E or finite displacement D, it can be switched to different fields 
        in different representations. If an electrical field is switched on the unit of the field 
        vector is V/Angstrom. If an dielectric displacement field is used the the unit 
        is e/Angstrom^2. For simulations with fluctuating cells it is recommended by SSV to used 
        reduced field variables. This is only possible together with a constand D simulation.
        In this case the field vector has units of charge e.
        
        See Nat Phys paper by SSV and the work by M. Sprijk

        :Parameters:
        """
        if reduced==False: assert const_D == True
        self.use_efield=True
        if self.is_master: self.pdlpio.rest_data.append("imgidx")
        # adjust units of field according to type of calc
        if const_D:
            self.set_calc_dipole(use_ref=False)
            self.set_atoms_moved()
            self.calc_energy_force()
            self.efield_const_D=True
            self.dlp.keyfld = 10
            ### put reference into the d part
            dip = self.get_dipole()
            V   = self.get_cell_volume()
            pol = dip/V
            ### calc reciprocal lattice vectors
            r_cell = self.get_reciprocal_cell()
            d_ref  = V*numpy.dot(r_cell, pol)
            print "pol:", pol
            self.pprint("Using constant electric displacement d [e]")
            # if not the reduced D is given we have to calculate it out of D
            if reduced == False:
                d = V*numpy.dot(r_cell,field_vect)
            else:
                d = field_vect
            self.pprint(field_vect)
            if use_ref:
                self.dlp_efld.prmfld[:3] = numpy.array(d_ref+d)
            else:
                self.dlp_efld.prmfld[:3] = numpy.array(d)
            self.set_atoms_moved()
        else:
            self.set_calc_dipole(use_ref,ref)
            self.efield_const_D=False
            self.dlp.keyfld = 1
            self.pprint("Using constant electric field E [V/Angstrom]")
            self.pprint(field_vect)
            # the constant 9648.53082 is the farday cosnt times 0.1 because DL_poly uses 10J/mol
            # as energy unit. we do not need to correct the length becasue the field is given per angstrom already
            self.dlp_efld.prmfld[:3] = numpy.array(field_vect)*9648.53082
            self.set_atoms_moved()
        return
        
    def change_efield(self, field_vect):
        """
        Use this if the field is already on and you just want to change the vector
        NOTE: you can not change the mode between const_D and const_E
        """
        assert self.use_efield==True, "Field needs to be switched on with set_efield() first"
        if self.efield_const_D:
            self.pprint("Using constant electric displacement D [e/Angstrom^2]")
            self.pprint(field_vect)
            self.dlp_efld.prmfld[:3] = numpy.array(field_vect)
            self.set_atoms_moved()
        else:
            self.pprint("Using constant electric field E [V/Angstrom]")
            self.pprint(field_vect)
            # the constant 9648.53082 is the farday cosnt times 0.1 because DL_poly uses 10J/mol
            # as energy unit. we do not need to correct the length becasue the field is given per angstrom already
            self.dlp_efld.prmfld[:3] = numpy.array(field_vect)*9648.53082
            self.set_atoms_moved()
        return

    def unset_efield(self):
        # switch off electric field
        self.pprint("Switching electric field OFF")
        dlp.keyfld = 0
        self.use_efield=False
        self.dlp_efld.prmfld[:] = 0.0
        self.dlp_efld.efield_energy = False
        return
        
    def set_calc_dipole(self, use_ref=True, ref = None):
        # switch on dipole/polariaztion calcualtion
        self.pprint("Switching on the computation of polarization/cell dipole")
        if not self.trackimg:
            self.track_images()
        self.dlp_efld.lcalc_dip = True
        self.dlp_efld.dip_ref[:] = 0.0
        self.set_atoms_moved()
        if use_ref:
            self.pprint("The initial reference dipole in e*A is:")
            if ref is None:                
                self.calc_energy()
                ref_dip = self.get_dipole()
                self.pprint(ref_dip)
                self.dlp_efld.dip_ref[:] = ref_dip
            else:
                self.pprint(ref)
                self.dlp_efld.dip_ref[:] = numpy.array(ref)
        return

    def get_dipole(self):
        if not self.dlp_efld.lcalc_dip:
            return numpy.zeros([3], "d")
        dipole = self.dlp_efld.dip
        return dipole
        
    def track_images(self, imgidx = None):
        if self.trackimg == False:
            self.dlp_conf.alloc_trackimg_arrays(self.idnode, self.nodes)
            self.trackimg = True
            self.pprint("Tracking of perioidc images: arrays allocated")
        self.dlp_conf.ltrackimg = True
        self.pprint("Tracking of perioidc images of atoms switched on")
        if imgidx is not None:
            self.dlp_conf.imgx = imgidx[:,0]
            self.dlp_conf.imgy = imgidx[:,1]
            self.dlp_conf.imgz = imgidx[:,2]
        return
        
    def track_images_off(self):
        self.dlp_conf.ltrackimg = False
        self.pprint("Tracking of perioidc images of atoms switched off")
        return


    def set_double_wall(self, split_index, z1, V1, z2, V2):
        self.pprint("Harmonic Double Wall Potential switched on")
        self.pprint("First set of atoms  (atom 1 to atom %d)" % (split_index-1))
        self.pprint("Bound to range between z1 %12.6f and z2 %12.6f" % (z1, z2))
        self.pprint("Second set of atoms (atom %d to atom %d)" % (split_index, self.natoms))
        self.pprint("Bound to the other part of the cell")
        self.pprint("Force constant: V1 %12.6f kcal/mol/A**2 , V2 %12.6f kcal/mol/A**2" % (V1, V2))
        if dlp.imcon>2:
            raise ValueError, "Double Wall Potential not for triclinic cells"
        if dlp.imcon==0:
            raise ValueError, "Double Wall Potential not for nonperidoic systems"
        cell = self.get_cell()
        zmax = cell[2,2]
        if (z1>z2) or (z1>zmax) or (z2>zmax) or (z1<0.0) or (z2<0.0):
            raise ValueError, "Inconsistent settings in wall positions!"
        dlp.keyfld = 8
        self.efld.prmfld[0] = z1
        self.dlp_efld.prmfld[1] = V1*self.dlp.engunit
        self.dlp_efld.prmfld[2] = z2
        self.dlp_efld.prmfld[3] = V2*self.dlp.engunit
        self.dlp_efld.split_index = split_index
        return

    def set_gaussian_wall(self, z1, V1, s1, z2, V2, s2):
        self.pprint("Gaussian Wall potential switched on")
        dlp.keyfld = 9
        self.dlp_efld.prmfld[0] = z1
        self.dlp_efld.prmfld[1] = V1*self.dlp.engunit
        self.dlp_efld.prmfld[2] = s2
        self.dlp_efld.prmfld[3] = z2
        self.dlp_efld.prmfld[4] = V2*self.dlp.engunit
        self.dlp_efld.prmfld[5] = s2
        return


    ######################################################################
    #                basic manipulation routines for changing
    #                atom positions/cell parameter etc
    ######################################################################

    def get_natoms(self):
        """ get the systems number of atoms """
        return self.natoms

    def get_tstep(self):
        """ get the timestep in ps """
        return self.dlp.tstep

    def get_elements(self):
        """ get a list of all elements"""
        return self.mol.elems

    def get_atomtypes(self, internal=False):
        """ get a list of the atomtypes """
        atypes = []
        for i in xrange(self.natoms):
            if internal:
                atypes.append(string.strip(self.dlp.get_atmnam(i+1)))
            else:
                atypes.append(self.mol.types[i])
        return atypes

    def get_xyz(self):
        """
        get the atomic positions as a numpy array
        for convenience we remap the linear coordinates in a [N,3] shape
        """
        xyz = numpy.empty([self.natoms, 3],"d")
        xyz[:,0] = self.dlp_conf.xxx
        xyz[:,1] = self.dlp_conf.yyy
        xyz[:,2] = self.dlp_conf.zzz
        return xyz

    def get_xyz_img(self):
        """
        get the atomic positions including the image info (wrt to when track imag was switched on)
        """
        if (not self.trackimg) or self.dlp.imcon>3:
            return self.get_xyz()
        if self.dlp.imcon == 3:
            # project into fractional coordinates
            xyz = self.get_xyz()
            cell = self.get_cell()
            inv_cell = numpy.linalg.inv(cell)
            frac_xyz = numpy.dot(xyz, inv_cell)
        #    print numpy.sum(abs(self.dlp_conf.imgx))
            frac_xyz[:,0] += self.dlp_conf.imgx
            frac_xyz[:,1] += self.dlp_conf.imgy
            frac_xyz[:,2] += self.dlp_conf.imgz
            xyz = numpy.dot(frac_xyz, cell)
            return xyz
        else:
            cell = self.get_cell().diagonal()
            xyz = numpy.empty([self.natoms, 3],"d")
            xyz[:,0] = self.dlp_conf.xxx + cell[0]*self.dlp_conf.imgx
            xyz[:,1] = self.dlp_conf.yyy + cell[1]*self.dlp_conf.imgy
            xyz[:,2] = self.dlp_conf.zzz + cell[2]*self.dlp_conf.imgz
        return xyz

    def get_imgidx(self):
        """
        Method to get the imgidx for constant E/D simulations
        """
        if not self.use_efield: return numpy.zeros([self.natoms,3],"i")
        imgidx = numpy.empty([self.natoms,3],"i")
        imgidx[:,0] = self.dlp_conf.imgx
        imgidx[:,1] = self.dlp_conf.imgy
        imgidx[:,2] = self.dlp_conf.imgz
        return imgidx


    def get_subset_xyz(self, aind):
        """ get atom positions for a subset as a numpy array
        :Parameters:

            - aind : list of atom indices (starting with 0, Python style)
        """
        xyz = numpy.empty([len(aind), 3], "d")
        xyz[:,0] = self.dlp_conf.xxx[aind]
        xyz[:,1] = self.dlp_conf.yyy[aind]
        xyz[:,2] = self.dlp_conf.zzz[aind]
        return xyz

    def set_xyz(self, new_xyz):
        """ set all atom positions from numpy array
        :Parameters:

            - new_xyz [natoms, 3] : numpy array with cartesian coordinates in Angstrom
        """
        self.dlp_conf.xxx[:] = new_xyz[:,0]
        self.dlp_conf.yyy[:] = new_xyz[:,1]
        self.dlp_conf.zzz[:] = new_xyz[:,2]
        self.set_atoms_moved()
        return

    def set_subset_xyz(self, new_xyz, aind):
        """ same as ``set_xyz`` for a subset
        :Parameters:

            - new_xyz [len(aind), 3] : numpy array with cartesian coordinates in Angstrom
            - aind : list of atom indices (starting with 0, Python style)

          """
        self.dlp_conf.xxx[aind] = new_xyz[:,0]
        self.dlp_conf.yyy[aind] = new_xyz[:,1]
        self.dlp_conf.zzz[aind] = new_xyz[:,2]
        self.set_atoms_moved()
        return

    def get_vel(self):
        """ returns velocities (internal units Angstrom per picosecond) as a numpy array """
        vel = numpy.empty([self.natoms, 3],"d")
        vel[:,0] = self.dlp_conf.vxx
        vel[:,1] = self.dlp_conf.vyy
        vel[:,2] = self.dlp_conf.vzz
        return vel

    def get_subset_vel(self, aind):
        """ returns velocities (internal units Angstrom per picosecond) as a numpy arrays """
        vel = numpy.empty([len(aind), 3],"d")
        vel[:,0] = self.dlp_conf.vxx[aind]
        vel[:,1] = self.dlp_conf.vyy[aind]
        vel[:,2] = self.dlp_conf.vzz[aind]
        return vel

    def set_vel(self, new_vel):
        """ set velocities from numpy array """
        self.dlp_conf.vxx[:] = new_vel[:,0]
        self.dlp_conf.vyy[:] = new_vel[:,1]
        self.dlp_conf.vzz[:] = new_vel[:,2]
        return

    def set_subset_vel(self, new_vel, aind):
        """ set subset velocites from numpy array (see set_subste_xyz) """
        self.dlp_conf.vxx[aind] = new_vel[:,0]
        self.dlp_conf.vyy[aind] = new_vel[:,1]
        self.dlp_conf.vzz[aind] = new_vel[:,2]
        return

    def set_atoms_moved(self):
        """
        tells the system that the energy/forces are invalid and need to be recomputed.
        most methods like :class:`pydlpoly.set_xyz` do this automatically
        So it is usually not necessary to call this from the user side
        """
        self.atoms_moved = True
        if self.qeq:
            self.qeq.set_recalc()
            if self.use_Jij: self.valid_Jij= False
        return

    def get_cell(self):
        """ get curernt cell vectors as a numpy array [3,3] """
        if (self.dlp.imcon == 0):
            self.pprint ("No cell: not a periodic system!")
            return
        cell = numpy.array(self.dlp_conf.cell,"d")
        return cell.reshape((3,3))

    def get_reciprocal_cell(self):
        """
        Return the reciprocal lattice multiplied by a factor of 2pi
        """
        cell = self.get_cell()
        return numpy.linalg.inv(cell).T

    def set_cell(self, cell, cell_only=False):
        """
        set cell from cell vectors

        :Parameters:

            - cell [3,3] : numpy array with cell vectors
            - cell_only : boolean if true only the cell is changed but the atom postions are untouched. \
                          by default the atoms are scaled according to the change of the cell
        """
        if not cell_only:
            abc = self.get_frac()
            # now recompute new xyz using the new cell
            new_xyz = numpy.dot(abc,cell)
            self.set_xyz(new_xyz)
        cell = cell.ravel()
        self.dlp_conf.cell[:] = cell[:]
        self.dlp.update_celprp()
        self.set_atoms_moved()
        return

    def get_frac(self):
        """
        get the current coordinates as fractional coordinates
        """
        # for this we need the current cell, invert it
        cell = self.get_cell()
        inv_cell = numpy.linalg.inv(cell)
        xyz = self.get_xyz()
        # compute fractional coordinates
        frac = numpy.dot(xyz,inv_cell)
        return frac

    def set_frac(self, frac):
        """
        set coordinates from fractional coords
        """
        cell = self.get_cell()
        new_xyz = numpy.dot(frac,cell)
        self.set_xyz(new_xyz)
        self.set_atoms_moved()
        return

    def get_cell_volume(self):
        """ get cell volume """
        if (self.dlp.imcon == 0):
            self.pprint ("No cell: not a periodic system!")
            return
#        not sure if the stpval entries (averaged?) still exist
#        return dlp_prop.stpval[18]
        if self.dlp.stpvol == 0.0:
            self.dlp.update_celprp()
            return self.dlp_conf.celprp[9]
        else:
            return self.dlp.stpvol

    def get_stress(self):
        """ get the stress tensor as a numpy array [3,3] in kcal/mol
        CAUTION: NOT the physical stress tensor, has to be divided by the Volume (use get_stress_tensor)
        """
        if (self.dlp.imcon == 0):
            self.pprint ("No stress: not a periodic system!")
            return
        stress = numpy.array(self.dlp_conf.stress,"d")
        return stress.reshape((3,3))/self.dlp.engunit

    def get_stress_tensor(self):
        """ get the stress tensor as a numpy array [3,3] in kcal/mol/A^3 """
        if self.dlp.imcon == 0:
            self.ppritn("No stress: not a periodic system")
            return
        else:
            return self.get_stress() / self.get_cell_volume()

    def get_temperature(self):
        """ get the (reference) system temperature """
        return self.dlp.temp

    def get_bcond(self):
        """ get the boundary conditions as an integer

        * 0 : non-periodic
        * 1 : cubic
        * 2 : orthormobic
        * 3 : triclinic
        * 6 : 2D periodic

        for higher numbers see ``DL_Poly`` documentation
        """
        return self.dlp.imcon

#RS WARNING this is commented out because with the recent changes the statistics are no longer calculated
#           for Michaels MetaMD we need to reimplement our own stattistics for the pressure tensor
    #def get_pressure_tensor(self):
        ## this returns the pressure tensor in dlpoly pressure unit prsunt
        ## as accumulated in the stpval and sumval arrays in the property_module
        ## the sumval is averaged over numacc steps
        #if (dlp.imcon == 0):
            #self.pprint ("No stress: not a periodic system!")
            #return
        #pres      = dlp_prop.stpval[27+dlp.ntpatm:27+9+dlp.ntpatm]
        #avrg_pres = dlp_prop.sumval[27+dlp.ntpatm:27+9+dlp.ntpatm]
        #pres.shape=(3,3)
        #avrg_pres.shape=(3,3)
        #return (pres, avrg_pres, dlp.numacc)

    def get_cellforce(self):
        cell    = self.get_cell()
        stress  = self.get_stress()
        # compute force from stress tensor
        cell_inv = numpy.linalg.inv(cell)
        cellforce = numpy.dot(cell_inv, stress)
        # at this point we could constrain the force to maintain the boundary conditions
        if self.dlp.imcon ==  2:
            # orthormobic: set off-diagonal force to zero
            cellforce *= numpy.eye(3)
        elif self.dlp.imcon == 6:
            # 2D periodic: remove off-diagonal terms to maintain orthorombic box and z component
            cellforce *= numpy.eye(3)
            cellforce[2,2] = 0.0
        if self.dlp.imcon == 1:
            # in the cubic case average diagonal terms
            avrgforce = cellforce.trace()/3.0
            cellforce = numpy.eye(3)*avrgforce
        return cellforce

    # part of the API for QEq and related fluctuating charge models

    def get_charges(self):
        """
        get the atomic charges as numpy array [natoms]
        """
        return numpy.array(self.dlp_conf.chge[:])

    def set_charges(self, q):
        """ set the atomic charges """
        self.dlp_conf.chge[:] = q[:]
        return

    def set_onsite_charge_energy(self, Ec_ii):
        # sets the contribution of the onsite charge terms (added to Ecoul during energy calc)
        # divide by nodes in order to compensate the global sum in parallel
        self.dlp_conf.enc_self = Ec_ii*self.dlp.engunit/float(self.nodes)
        return

    def get_sigmas(self):
        return numpy.array(self.dlp_conf.qsigma[:])

    def set_sigmas(self, sigmas):
        self.dlp_conf.qsigma[:] = sigmas[:]
        return

    def get_masses(self):
        """ get atomic masses """
        return numpy.array(self.dlp_conf.weight[:])

    def get_subset_masses(self, aind):
        """ get masses for atoms indexed in list aind """
        return numpy.array(self.dlp_conf.weight[aind])
        
    def unwrap_molecules(self):
        xyz = self.get_xyz()
        cell = self.get_cell()
        celldiag = numpy.diagonal(cell)
        for i in range(self.mol.nmols):
            mol_xyz = xyz[self.mol.mols[i]]
            dxyz = mol_xyz-mol_xyz[0]
            dxyz -= celldiag*numpy.around(dxyz/celldiag)
            dxyz += mol_xyz[0]
            xyz[self.mol.mols[i]] = dxyz[:]
        self.set_xyz(xyz)
        return
        
    def unwrap_fragments(self, frags):
        xyz = self.get_xyz()
        cell = self.get_cell()
        celldiag = numpy.diagonal(cell)
        for f in frags:
            frag_xyz = xyz[f]
            dxyz = frag_xyz-frag_xyz[0]
            dxyz -= celldiag*numpy.around(dxyz/celldiag)
            dxyz += frag_xyz[0]
            xyz[f] = dxyz[:]
        self.set_xyz(xyz)
        return

    def get_loc_natoms(self):
        # this is clumsy but safe
        self.loc_natoms=0
        for i in xrange(self.idnode,self.natoms,self.nodes):self.loc_natoms+=1
        return self.loc_natoms

    def get_loc_listsize(self):
        if not self.loc_natoms: self.get_loc_natoms()
        self.loc_nlist = numpy.sum(self.dlp_conf.lentry[:self.loc_natoms])
        self.loc_exlist= numpy.sum(self.dlp_excl.nexatm[:self.loc_natoms])
        # print self.idnode, self.loc_natoms, self.loc_nlist, self.loc_exlist
        return self.loc_nlist, self.loc_exlist

    def inquire_core_charge_use(self):
        return self.dlp_conf.use_core_charge

    def get_core_charges(self):
        return numpy.array(self.dlp_conf.core_chge[:])

    def get_data_funcs(self):
        """ returns a directory of names and corresponding functions
            the functions return pure double precision numpy arrays,
            which can be stored by pdlpio as restart or trajectory data"""
        data_funcs = {"xyz"      : self.get_xyz,\
                      "vel"      : self.get_vel,\
                      "charges"  : self.get_charges,\
                      "forces"   : self.get_force,\
                      "xyz_img"  : self.get_xyz_img,\
                      }
        if self.dlp.imcon > 0:
            data_funcs["cell"]     = self.get_cell
            data_funcs["cellfrcs"] = self.get_cellforce
            data_funcs["stress"]   = self.get_stress_tensor
            data_funcs["imgidx"]   = self.get_imgidx
        data_funcs["dipole"]   = self.get_dipole
        for sp in self.scal_prop:
            data_funcs[sp] = lambda sp=sp: self.get_scal_prop(sp)
        return data_funcs

    def get_conn(self):
        return self.mol.cnct

    ######################################################################
    #                thermostating
    #
    ######################################################################

    def use_plumed_metamd(self, plumedinput):
        """ this will switch on the use of plumed metadynamic in pydlpoly
        NOTE: this is always active ... do not use in optimizations
        """
        if not plumed_available:
            self.pprint("PyPlumed Module is not available!!")
            self.pprint("ABORTING")
            return
        if self.plumed_init:
            self.pprint("Plumed is already initialized !!")
            self.pprint("ABORTING")
            return
        self.pprint("\nINIT PYPLUMED\n")
        # init pyplumed directly from dlpoly datat structures
        _pyplumed.init_metadyn_(self.dlp.tstep, self.dlp_conf.weight, self.dlp_conf.chge, self.dlp.imcon, self.dlp.engunit, plumedinput, 0)
        self.plumed_init = True
        self.plumed_step = 0
        self.meta_eng = 0.0
        return

    def plumed_energy(self):
        self.plumed_step += 1
        """ this method is called in the energy function every step if metamd is init and on """
        metaenergy = _pyplumed.meta_force_calculation_(self.dlp_conf.cell, self.plumed_step,\
                                     self.dlp_conf.xxx, self.dlp_conf.yyy, self.dlp_conf.zzz,\
                                     self.dlp_conf.fxx, self.dlp_conf.fyy, self.dlp_conf.fzz)
        return metaenergy/self.dlp.engunit

    def metamd_on(self):
        self.plumed=True
        return

    def metamd_off(self):
        self.plumed=False
        self.plumed_step = 0
        return

    ######################################################################
    #                thermostating
    #
    ######################################################################

    def set_thermostat_sigma(self, degfree):
        """ resets the sigma value used in the thermostat to determine the reference kin energy"""
        self.sigma = 0.5*degfree*boltz*dlp.temp
        self.dlp.sigma  = self.sigma
        return

    def get_degfree(self):
        """ get degrees of dreedom from Fortran"""
        return self.dlp.degfre

    def set_degfree(self, dof):
        """ set degrees of freedom and sigma (for thermostat) in fortran """
        self.dlp.degfre = float(dof)
        self.sigma = 0.5*self.dlp.degfre*boltz*self.dlp.temp
        self.dlp.sigma  = self.sigma
        return

    ######################################################################
    #                basic MD routines
    #
    ######################################################################

    def MD_init(self, stage, T=None, p=None, startup=False, ensemble="nve", thermo=None, relax=None, \
                             traj=None, rnstep=100, tnstep=100, grlmb=False, grlmb_N=0):
        """
        Sets up everything for a new MD stage (you MUST provide a new stage name here)

        :Parameters:
            - stage (str)    : name of the new stage to be used in restart
            - T (float)      : temperature to be used for the thermostat (and for the intial veldist) [K]
            - p (float)      : pressure (in dlpoly units [kAtm]??check this??)
            - startup (bool) : if True use a Maxwell-Boltzmann distribution for the velocities given by T
            - ensemble (str) : type of ensemble (nve, nvt, npt, nst)
            - thermo (str)   : "ber" or "hoover" (includes barostat .. dlpoly manual)
            - relax (tuple or list of floats) : relaxation times (taut, taup) in [ps] , taup only for npt, nst
            - traj (list of str or tuple) : list of keys to keep in the trajectory file (if tuple: (key, nstep))
            - rnstep (int)   : number of steps after which restart info is written
            - tnstep (int or list of int)   : number of steps after which traj info is written (if list it needs to be of the same length as traj)
            - grlmb (`False` or "added" or molecule name or list of :class:`pdlpmol.pdlpmol` objects) : if this is not `False` the
                             self.added_mols will be switched 1 during the MD
            - grlmb_N (int)  : if grlmb is "added" or a molecule_name etc. and grlmb_N > 0 (not the default) then only the first
                               grlmb_N molecules are grown in and the others are kept at lambda=0
        """
        # check params
        # tnstep is a list ... then it must be of the same length as traj
        if traj:
            if type(tnstep) == types.ListType:
                assert len(tnstep) == len(traj)
                tnstep_list = tnstep
            else:
                # tnstep is the default value ...
                tnstep_list=len(traj)*[tnstep]
                for t in xrange(len(traj)):
                    if type(traj[t]) == types.TupleType:
                        tnstep_list[t] = traj[t][1]                    
                        traj[t] = traj[t][0]
        else:
            tnstep_list = []
        # start up MD stage
        self.md_counter = 0
        self.dlp.nstep=0
        self.dlp.numacc=0
        # set the ensemble type
        self.md_ens = ensemble.lower()
        if enstype.has_key(self.md_ens):
            keyens = enstype[self.md_ens]
        else:
            raise ValueError, "Ensemble %s is unknown" % ensemble
        # set theromstat settings
        self.md_thermo = thermo
        self.md_relax  = relax
        if self.md_ens != "nve":
            if thermo == "ber":
                pass
            elif thermo == "hoover":
                keyens += 1
            else:
                raise ValueError, "Unknown thermostat"
        self.dlp.keyens = keyens
        if ensemble[2] == "t":
            if not thermo :
                raise ValueError, "Thermostat must be set with thermo= "
            if not relax :
                raise ValueError, "You have to set thermostat/barostat relaxation params with relax="
            self.dlp.taut = relax[0]
            if ensemble[1] == "p" or ensemble[1] == "s":
                self.dlp.taup = relax[1]
        # set the temperature
        self.md_T = T
        if T:
            self.dlp.temp = T
        if startup:
            self.dlp.keyres = 0
        else:
            # set dlpoly restart flag to 1 = keep vleocities (this implies to read them from restart!!)
            self.dlp.keyres = 1
        # set pressure
        self.md_p = p
        if p:
            # we need to convert here from katm to dlpoly internal units
            self.dlp.press = p/prsunt
        # now call the startup procedure
        # NOTE: in the original dlpoly this is done only once but there also no multiple "stages" where possible
        #       the major work there is to zero some counters and statisitcs collectors.
        #       so it should be safe to call this multiple times but this needs to be tested!!!
        # NOTE 2: in init_dyn the degfre are determined this should be done on python level
        #         and be removed in fortran. this allows to change the number of degfre (this is a float)
        self.dlp.init_dyn(skip_REVOLD=True)
        # now make a new stage in the restart file and printout
        # currently only the master node has a working pdlpio object -> make pdlpio paralle and change this
        self.md_stage = stage
        self.md_rnstep = rnstep
        self.md_tnstep = tnstep_list
        self.md_traj   = traj
        if self.is_master:
            if self.pdlpio:
                if self.pdlpio.has_stage(stage):
                    raise IOError, "The PDLP restart file already has a stage of this name! Use recover"
                self.pdlpio.add_stage(self.md_stage, self.md_rnstep)
                if self.md_traj:
                    self.pdlpio.add_traj(self.md_tnstep, self.md_traj, tstep=self.get_tstep())
        self.pprint("\n\n>>> MD INIT >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        self.pprint("Starting up MD (setting counters back to zero)")
        self.pprint("   STAGE       %s" % stage)
        self.pprint("   ENSEMBLE    %s" % ensemble)
        if startup:
            self.pprint ("   STARTUP of atom velocities to Maxwell-Boltzmann distribution")
        if self.md_ens != "nve":
            out_thermo= "   THERMO      %s"  % self.md_thermo
            out_thermo+="        ("+(len(relax)*"%10.5f "%tuple(relax))+") [ps]"
            self.pprint(out_thermo)
            if not self.md_T: raise ValueError, "Ensemble requires a temperature to be set"
            self.pprint("   TEMP   %10.1f K" % self.md_T)
            if self.md_ens == "npt" or self.md_ens == "nst":
                if not self.md_p: raise ValueError, "Ensemble requires a pressure to be set"
                self.pprint("   PRES   %10.5f katm" % self.md_p)
        self.pprint("   REST steps  %d" % self.md_rnstep)
        if self.md_traj:
            if type(self.md_tnstep) == types.IntType:
                self.md_tnstep = [self.md_tnstep]
            self.pprint("   TRAJ writing   "+(len(self.md_traj)*"%8s " % tuple(self.md_traj)))
            self.pprint("        steps     "+(len(self.md_tnstep)*"%8d " % tuple(self.md_tnstep)))
        # LMABDA stuff
        if grlmb:
            self.pprint("   GROW LAMBDA: lamda will be set to zero and then switched to one linearly")
            if grlmb=="added":
                self.md_grlmb = self.added_mols
                self.pprint("   switching on added molecules")
            elif type(grlmb) == types.StringType:
                self.md_grlmb = self.molecules.get_list_from_name(grlmb)
                self.pprint("   switching on molecules with name %s" % grlmb)
            elif type(grlmb) == types.ListType:
                self.md_grlmb = grlmb
                self.pprint("   switching on molecules %s" % str(grlmb))
            else:
                raise ValueError, "Unknown option for grlmb"
            if self.use_split_vdw:
                for m in self.md_grlmb: m.set_lambda(0.0, 0.0, 0.0)
            else:
                for m in self.md_grlmb: m.set_lambda(0.0, 0.0)
            # check if only a subset should be grown in and reduce md_grlmb accordingly (all other lambdas have been set to zero
            # and will stay there
            if grlmb_N > 0:
                assert grlmb_N <= len(self.md_grlmb), "grlmb_N is larger then the available molecules"
                self.pprint("   growing in only the first %d molecules is requested" % grlmb_N)
                self.md_grlmb = self.md_grlmb[:grlmb_N]
        return


    def MD_run(self,nsteps, do_every_step=None, printout=100):
        self.mode="MD"
        self.pprint("performing %d MD steps" % nsteps)
        # NOTE: here sigma is calculated from temp and degfre
        self.dlp.dyns_start()
        label = "#%-5s" % (self.md_stage.upper()[:5])
        self.report_energy_header(marker=("            "))
        self.calc_energy()
        for i in xrange(nsteps):
            if ((i%printout)==0):
                self.report_energies(marker=("%s %8d "%(label,i)))
            self.timer.switch_on("MD prop fh")
            self.dlp.dyns_step_fh()
            for s in self.extra_systems: s.vv_fh()
            self.set_atoms_moved()
            self.timer.switch_off()
            self.calc_energy()
            self.timer.switch_on("MD prop sh")
            self.dlp.dyns_step_sh()
            for s in self.extra_systems: s.vv_sh()
            # do restart only on master (pdlpio not None)
            if self.pdlpio:
                self.timer.switch_to("MD io")
                self.pdlpio()
            # do this every step
            self.timer.switch_to("MD misc")
            if do_every_step: do_every_step(self, i, nsteps)
            # grow lambda option
            if self.md_grlmb:
                lmb = float(i+1)/float(nsteps)
                lmb *= lmb*lmb
                if self.use_split_vdw:
                    for m in self.md_grlmb: m.set_lambda(lmb, lmb, lmb)
                else:
                    for m in self.md_grlmb: m.set_lambda(lmb, lmb)
                if ((i%printout)==0):
                    self.pprint("         -> grow lambda lmb=%10.5f" % lmb)
            self.timer.switch_off()
        self.report_energies(marker=("%s %5d "%(label,i+1)))
        if self.pdlpio: self.pdlpio(force_wrest = True)
        if self.md_grlmb:
            self.pprint("         -> grow lambda -- switching all lambdas back to 1.0")
            if self.use_split_vdw:
                for m in self.md_grlmb: m.set_lambda(1.0, 1.0, 1.0)
            else:
                for m in self.md_grlmb: m.set_lambda(1.0, 1.0)
            # now switch md_grlmb to None in order to avoid any problems with lambda in consecutive MD steps
            self.md_grlmb = None
        # NOTE: this just calls -> result -> revive (final printout and write to REVOLD)
        # dlp.dyns_finish()
        return

    ######################################################################
    #                basic optimization routines
    ######################################################################


    def MIN_sd(self,thresh,maxiter=1000):
        self.mode ="MIN"
        step = 0.000002
        stop = False
        i = 0
        self.pprint("\n\nSTEEPEST DESCENT for %d steps (thresh=%10.5f)" % (maxiter, thresh))
        while not stop:
            i += 1
            coords = self.get_xyz()
            energy, force  = self.calc_energy_force()
            rms = numpy.sqrt(numpy.sum(force*force)/(3*self.natoms))
            if rms <= thresh : stop = True
            coords += step*force
            self.set_xyz(coords)
            if self.pdlpio: self.pdlpio()
            self.pprint ("%3d  %12.6f %12.6f " % (i, energy, rms))
            if i>maxiter:
                stop=True
                self.pprint("Maximum number of iterations reached! Not Converged")
        self.report_energies()
        return

    def MIN_dlpoly_cg(self,thresh,maxiter=1000):
        # use the dl_poly optimizer (rather poor CG optimizer with lot of problems in the inital phase)
        # we switch the loptim flag to true
        # this means the coordinates are uptdated any time the energy is called
        # therfore we need to set the self.atoms_moved flag to True every step
        selfmode="MIN"
        self.dlp.loptim = True
        stop = False
        i = 0
        self.pprint("\n\nDL_Poly Native CG Optimizer for %d steps (thresh=%10.5f)" % (maxiter, thresh))
        while not stop:
            i += 1
            self.set_atoms_moved()
            energy, force  = self.calc_energy_force()
            rms = numpy.sqrt(numpy.sum(force*force)/(3*self.natoms))
            if rms <= thresh : stop = True
            if self.pdlpio: self.pdlpio()
            self.pprint ("%3d  %12.6f %12.6f " % (i, energy, rms))
            if i>maxiter:
                stop=True
                self.pprint("Maximum number of iterations reached! Not Converged")
        self.report_energies()
        self.dlp.loptim = False
        return

    def optgeom_callback(self, x, g):
        """ callback for pure geometry optimizations
            both x and g are flat arrays
            energy and gradient are in engunit (kcal/mol)
        """
        self.dlp_conf.xxx[:] = x[0*self.natoms:1*self.natoms]
        self.dlp_conf.yyy[:] = x[1*self.natoms:2*self.natoms]
        self.dlp_conf.zzz[:] = x[2*self.natoms:3*self.natoms]
        self.set_atoms_moved()
        # test for neigbor list update
        if self.count_MIN_lbfgs%self.nlst_rebuild == 0:
            self.dlp_nlst.nlst_rebuild = True
        e = self.calc_energy()
        #self.pprint("DEBUG opt callback e %12.6f count %d rebuild %d" % (e, self.count_MIN_lbfgs, dlp_nlst.nlst_rebuild))
        self.dlp_nlst.nlst_rebuild = False
        g[0*self.natoms:1*self.natoms] = self.dlp_conf.fxx[:]
        g[1*self.natoms:2*self.natoms] = self.dlp_conf.fyy[:]
        g[2*self.natoms:3*self.natoms] = self.dlp_conf.fzz[:]
        g *= -1.0/self.dlp.engunit
        self.count_MIN_lbfgs += 1
        # write restart
        if self.pdlpio: self.pdlpio()
        return e


    def MIN_lbfgs(self, thresh, maxiter=None, m=10, nlst_rebuild=10, verbose = True):
        """
        Minimize the systems energy by a L-BFGS optimizer (no lattice optimization)

        :Parameters:
             - thresh:       (float) Threshold of RMS Gradient in kcal/mol/A
             - maxiter:      (int)   Maximum Number of iterations
             - m:            (int)   Number of previous gradients used in L-BFGS
             - nlst_rebuild: (int)   Number of steps until a rebuild of neigbor list is enforced
        """
        self.mode="MIN"
        self.pprint ("\n\n")
        self.pprint(self.mark_line)
        self.pprint ("LBFGS geometry optimization using a history of %d" % m)
        self.pprint ("    convergence threshold: %10.5f"%thresh)
        # generate an optimizer instance (note: threshold is handeled)
        # NEW interface for MExt: first argument is pydlpoly instance to get access to lbfgs_module
        lbfgs_opt = lbfgs.lbfgs(self, self.natoms*3, m, self.optgeom_callback, thresh, verbose, 
                mpi_comm = self.local_comm, out=self.out)
        # make a counter
        self.count_MIN_lbfgs = 0
        # switch off nlist testing
        self.dlp_nlst.nlst_notest = True
        # keep rebuild interval
        self.nlst_rebuild = nlst_rebuild
        # generate arrays
        x = numpy.empty([3*self.natoms], "d")
        g = numpy.empty([3*self.natoms], "d")
        # set initial values
        x[0*self.natoms:1*self.natoms] = self.dlp_conf.xxx[:]
        x[1*self.natoms:2*self.natoms] = self.dlp_conf.yyy[:]
        x[2*self.natoms:3*self.natoms] = self.dlp_conf.zzz[:]
        # optimize
        lbfgs_opt(x, g, maxiter=maxiter)
        # set final positions
        self.dlp_conf.xxx[:] = x[0*self.natoms:1*self.natoms]
        self.dlp_conf.yyy[:] = x[1*self.natoms:2*self.natoms]
        self.dlp_conf.zzz[:] = x[2*self.natoms:3*self.natoms]
        self.set_atoms_moved()
        # recalc final energy here
        energy, force  = self.calc_energy_force()
        #rms_force = numpy.sqrt(numpy.sum(force*force)/(3*self.natoms))
        #self.pprint ("  final energy %12.6f gradient norm %12.6f" % (energy, rms_force))
        if self.pdlpio : self.pdlpio(force_wrest=True)
        self.report_energies(full=True)
        self.pprint ("LBFGS done")
        self.pprint(self.mark_line)
        # switch nlist testing back on
        self.dlp_nlst.nlst_notest = False
        return

    def LATMIN_sd(self,threshlat, thresh, lat_maxiter= 100, maxiter=5000, fact = 2.0e-3, maxstep = 3.0):
        """
        Lattice and Geometry optimization (uses MIN_lbfgs for geom opt and steepest descent in lattice parameters)

        :Parameters:
            - threshlat (float)  : Threshold in RMS force on the lattice parameters
            - thresh (float)     : Threshold in RMS force on geom opt (passed on to :class:`pydlpoly.MIN_lbfgs`)
            - lat_maxiter (int)  : Number of Lattice optimization steepest descent steps
            - fact (float)       : Steepest descent prefactor (fact x gradient = stepsize)
            - maxstep (float)    : Maximum stepsize (step is reduced if larger then this value)

        """
        self.mode = "LATMIN"
        self.pprint ("\n\nLattice Minimization: using steepest descent for %d steps (threshlat=%10.5f, thresh=%10.5f)" % (lat_maxiter, threshlat, thresh))
        self.pprint ("                      the geometry is relaxed with LBFGS at each step for a mximum of %d steps" % maxiter)
        self.pprint ("Initial Optimization ")
        self.MIN_lbfgs(thresh, maxiter=maxiter)
        # to zero the stress tensor
        self.dlp.nstep = 1
        self.set_atoms_moved()
        oldenergy = self.calc_energy()
        cell = self.get_cell()
        self.pprint ("Initial cellvectors:\n%s" % numpy.array2string(cell,precision=4,suppress_small=True))
        cellforce = self.get_cellforce()
        self.pprint ("Initial cellforce:\n%s" % numpy.array2string(cellforce,precision=4,suppress_small=True))
        stop = False
        latiter = 1
        while not stop:
            self.pprint ("Lattice optimization step %d" % latiter)
            step = fact*cellforce
            steplength = numpy.sqrt(numpy.sum(step*step))
            self.pprint("Unconstrained step length: %10.5f Angstrom" % steplength)
            if steplength > maxstep:
                self.pprint("Constraining to a maximum steplength of %10.5f" % maxstep)
                step *= maxstep/steplength
            new_cell = cell + step
            self.pprint ("New cell:\n%s" % numpy.array2string(new_cell,precision=4,suppress_small=True))
            self.set_cell(new_cell)
            self.MIN_lbfgs(thresh, maxiter=maxiter)
            # to zero the stress tensor
            self.dlp.nstep = 1
            self.set_atoms_moved()
            energy = self.calc_energy()
            if energy > oldenergy:
                self.pprint("WARNING: ENERGY SEEMS TO RISE!!!!!!")
            oldenergy = energy
            cell = self.get_cell()
            #self.pprint ("Current cellvectors:\n%s" % str(cell))
            cellforce = self.get_cellforce()
            self.pprint ("Current cellforce:\n%s" % numpy.array2string(cellforce,precision=4,suppress_small=True))
            rms_cellforce = numpy.sqrt(numpy.sum(cellforce*cellforce)/9.0)
            self.pprint ("Current rms cellforce: %12.6f" % rms_cellforce)
            latiter += 1
            if latiter >= lat_maxiter: stop = True
            if rms_cellforce < threshlat: stop = True
        self.pprint ("SD minimizer done")
        return


    ######################################################################
    #             Methods for ASE added by JPD
    ######################################################################

    def create_ase_atoms_object(self):
        from ase import Atoms
        a_elements = []
        for i in self.get_elements():
            e = string.split(i)[0]
            # HACK : make this general !! RS
            if e == 'CU':
                e ='Cu'
            #a_elements.append(string.split(i)[0])
            a_elements.append(e)
        a = Atoms(a_elements)
        a.set_positions(self.get_xyz())
        a.set_masses(self.get_masses())
        a.set_charges(self.get_charges())
        a.set_connectivity(self.mol.cnct)
        a.set_atomtypes(self.get_atomtypes())
        if (self.dlp.imcon != 0):
            a.set_cell(self.get_cell())
        return a

    def read_tinker_xyz(self, fname):
        xyz = []
        f = open(fname, "r")
        string.split(f.readline())
        for i in xrange(self.get_natoms()):
            lbuffer = string.split(f.readline())
            xyz.append(map(string.atof, lbuffer[2:5]))
        self.set_xyz(numpy.array(xyz))
        return

    ######################################################################
    #             Methods for TAMKIN added by JPD
    ######################################################################


    # not working due to commonblocks --> pickle to prep container
    def create_TAMKIN_molecule(self, hessian):
        from tamkin.data import Molecule
        import molmod
        positions = self.get_xyz() * ang2bohr
        atomicnumbers = []
        masses = self.get_masses() * g2au
        energy, force = self.calc_energy_force()
        energy *= kcal2au
        gradient = -kcal2au/ang2bohr * force
        hessian *= kcal2au/(ang2bohr**2)
        if (self.dlp.imcon != 0):
            cell = self.get_cell()*ang2bohr
            is_periodic = True
        else:
            is_periodic = False
        for i in range(self.get_natoms()):
            atomicnumbers.append(elements.labels.index(string.split(string.lower(self.get_elements()[i]))[0]))
        return Molecule(atomicnumbers, positions, masses, energy,
                gradient, hessian, 1,periodic=is_periodic,
                unit_cell = molmod.unit_cells.UnitCell(cell))



    ######################################################################
    #             Virtual Atoms
    #
    #  virtual atoms are of elemnt XX and have no mass, charge or vdw interactions
    #  (you need to make sure that this is the case in the key file)
    #
    #  they are positioned on the COM before the energy is computed and
    #  the forces are then distributed back to the defining atoms
    #
    ######################################################################

    def init_virt_atoms(self):
        self.virt_atom_scales = []
        mass = self.get_masses()
        for i, j in enumerate(self.virtual_atoms):
            # the scale is the atommass divided by the total mass of the fragment
            atoms = self.virt_atom_defs[i]
            m = mass[atoms]
            scale = m/m.sum()
            self.virt_atom_scales.append(scale)
        self.virt_atom_set_pos()
        return

    def virt_atom_set_pos(self):
        xyz = self.get_xyz()
        # precompute cell related things .. if bcond == 3 work in fractional coords
        cell = self.get_cell()
        if self.dlp.imcon == 3:
            icell = numpy.linalg.inv(cell)
        elif self.dlp.imcon > 0:
            celldiag = cell.diagonal()
            half_celldiag = 0.5*celldiag
        else:
            pass
        # do it parallel accumulate com an broadcast before
        nvirt = len(self.virtual_atoms)
        com = numpy.zeros([nvirt, 3], dtype="float64")
        for i in xrange(self.idnode,nvirt,self.nodes):
            atoms = self.virt_atom_defs[i]
            scale = self.virt_atom_scales[i]
            pos = xyz[atoms]
            # check for PBC
            if self.dlp.imcon == 3:
                fpos = numpy.dot(pos, icell)
                dfpos = fpos[1:]-fpos[0]
                fpos[1:] += numpy.where(numpy.less_equal(dfpos, -0.5), 1.0, 0.0)
                fpos[1:] -= numpy.where(numpy.greater   (dfpos,  0.5), 1.0, 0.0)
                pos = numpy.dot(fpos, cell)
            elif self.dlp.imcon > 0:
                dpos = pos[1:]-pos[0]
                pos[1:] += numpy.where(numpy.less_equal(dpos, -half_celldiag), 1.0, 0.0)*celldiag
                pos[1:] -= numpy.where(numpy.greater   (dpos,  half_celldiag), 1.0, 0.0)*celldiag
            else:
                pass
            com[i] = numpy.sum(pos*scale[:,None], axis=0)
        # "sum" up . essintially a broadcast
        com = self.gdsum(com)
        xyz[self.virtual_atoms] = com
        self.set_xyz(xyz)
        return

    def move_virtual(self,virtual,step):
        if virtual in self.virtual_atoms:
            atoms = self.virt_atom_defs[self.virtual_atoms.index(virtual)]
            xyz = self.get_xyz()
            xyz[virtual,:] += step
            for i in atoms:
                xyz[i,:] += step
            self.set_xyz(xyz)
        else:
            print '%s is not a virtual site' % virtual
            raise IOError
        return

    def virt_atom_distribute_force(self):
        """ for larger systems it might be better to directly access the fortran arrays fxx/fyy/fzz
            or to move this to fortran into a module """
        add_fxyz = numpy.zeros([self.natoms, 3], dtype="float64")
        fxyz = self.get_force()
        nvirt = len(self.virtual_atoms)
        # print numpy.sum(fxyz, axis=0)
        for i in xrange(self.idnode, nvirt, self.nodes):
            j = self.virtual_atoms[i]
            atoms = self.virt_atom_defs[i]
            scale = self.virt_atom_scales[i]
            vforce = fxyz[j]
            add_fxyz[atoms] += vforce[None,:]*scale[:,None]
            add_fxyz[j] = -vforce
            # zero the velocity of the virtual atom just to be safe
            #self.dlp_conf.vxx[j] = 0.0
            #self.dlp_conf.vyy[j] = 0.0
            #self.dlp_conf.vzz[j] = 0.0
        add_fxyz = self.gdsum(add_fxyz)
        self.add_force(add_fxyz)
        ### for frozen virtual sites ###
        if True in self.frozen_virtuals:
            for i in range(nvirt):
                nettoforce = numpy.zeros([3], dtype="float64")
                fxyz = self.get_force()
                add_fxyz = numpy.zeros([self.natoms, 3], dtype="float64")
                if self.frozen_virtuals[i] == True:
                    j = self.virtual_atoms[i]
                    atoms = self.virt_atom_defs[i]
                    scale = self.virt_atom_scales[i]
                    #print atoms
                    #print scale
                    for k in range(len(atoms)):
                        nettoforce += fxyz[atoms[k]]
                    for k in range(len(atoms)):
                        add_fxyz[atoms[k],:] = -nettoforce * scale[k]
                    #print self.get_force()
                    self.add_force(add_fxyz)
                    #print add_fxyz
            #print self.get_force()
        # DEBUG
        #print numpy.sum(self.get_force(),axis=0)
        return


    def project_force_on_virt(self, r_xyz, r_f):
        #print '#########################'
        nvirt = len(self.virtual_atoms)
        virt_xyz = numpy.zeros((nvirt,3),dtype = 'float64')
        virt_f = numpy.zeros((nvirt,3),dtype = 'float64')
        for i in range(nvirt):
            j = self.virtual_atoms[i]
            atoms = self.virt_atom_defs[i]
            virt_xyz[i,:] = r_xyz[j,:]
            for k in range(len(atoms)):
                #print r_f[atoms[k]]
                virt_f[i,:] += r_f[atoms[k]]
        #print r_f
        #print virt_f
        #print '#########################'
        return virt_xyz, virt_f


    def calc_nettoforce(self, virt):
        if virt in self.virtual_atoms:
            nettoforce = numpy.zeros([3], dtype="float64")
            fxyz = self.get_force()
            atoms = self.virt_atom_defs[self.virtual_atoms.index(virt)]
            for i in range(len(atoms)):
                nettoforce += fxyz[atoms[i],:]
        else:
            print "%s is not a virtual atom"
            raise IOError
        return nettoforce

    ######################################################################
    #     Utility functions ... in particular printout etc.
    ######################################################################

    def write_xyz(self, fname, add_info=None):
        """ writing out a simple xyz file for DEBUG purposes
        add_info: add charges or atomtypes for molden. 'Label -> atom+charges'."""
        if self.is_master:
            xyz = self.get_xyz()
            elms = self.get_elements()
            fxyz = open(self.rundir+"/"+fname,"w")
            fxyz.write("%d\n" % self.natoms)
            fxyz.write("PyDLP xyz file\n")
            if add_info == "atomtypes":
                info = self.get_atomtypes()
                for i in xrange(self.natoms):
                    infoi = float(''.join(c for c in info[i] if c.isdigit()))
                    fxyz.write("%s %12.6f %12.6f %12.6f %12.6f\n" % (elms[i], xyz[i,0], xyz[i,1], xyz[i,2], infoi))
            elif add_info == "charges":
                info = self.get_charges()
                for i in xrange(self.natoms):
                    infoi = info[i]
                    fxyz.write("%s %12.6f %12.6f %12.6f %12.6f\n" % (elms[i], xyz[i,0], xyz[i,1], xyz[i,2], infoi))
            elif type(add_info) == str:
                raise ValueError("add_info cannot be a different string than 'charges' or 'atomtypes'")
            elif add_info is None:
                for i in xrange(self.natoms):
                    fxyz.write("%s %12.6f %12.6f %12.6f\n" % (elms[i], xyz[i,0], xyz[i,1], xyz[i,2]))
            else:
                assert len(add_info) == self.natoms[()]
                for i in xrange(self.natoms):
                    infoi = add_info[i]
                    fxyz.write("%s %12.6f %12.6f %12.6f %12.6f\n" % (elms[i], xyz[i,0], xyz[i,1], xyz[i,2], infoi))
            fxyz.close()
        return

    def get_mollist(self, molname):
        moltype = self.mol.molnames.index(molname)
        return [i for i, n in enumerate(self.mol.moltypes) if n == moltype]

    def write_extxyz(self, fname, mols=None, nmols=None, zupfrac=None, rm_spc=[], fullcell=True):
        """ Writing out a simple xyz file

        :param fname: file name
        :param mols: molecules to write
        :param nmols: number of molecules to write sorted in z-ascending order
        :param zupfrac: fraction of z-length to choose the molecules
        :param rm_spc: do not write this species
        :param fullcell: write full cell
        :returns:
        :rtype:

        """
        xyz = self.get_xyz()
        cell = self.get_cell()

        def extxyz_head():
            if cell.any():
                if fullcell:
                    fields = ['Lattice=', 'Properties=', 'pbc=']
                    fields[0] += '"' + ' '.join(' '.join(str(cell) for cell in row) for row in cell) + '"'
                else:
                    fields = ['Cellpars=', 'Properties=', 'pbc=']
                    cellparams = unit_cell.abc_from_vectors(cell)
                    fields[0] += '"' + ' '.join(str(param) for param in cellparams) + '"'
                    rotcell = unit_cell.vectors_from_abc(cellparams)
                    # compute fractional coords in old cell
                    inv_cell = numpy.linalg.inv(cell)
                    abc = numpy.dot(xyz, inv_cell)
                    # compute cartesian coordinates in new cellvectors from fractional coords
                    xyz = numpy.dot(abc, rotcell)
                fields[1] += 'species:S:1:pos:R:3'
                bc = self.get_bcond()
                if bc == 0:
                    bc = 'F F F'
                elif bc == 6:
                    bc = 'T T F'
                else:
                    bc = 'T T T'
                fields[2] += '"' + bc + '"'
                title_line = " ".join(fields)
            else:
                title_line = 'PyDLP xyz file'
            return title_line

        if self.is_master:
            if mols:
                atomlist = []
                for m in mols:
                    mollist = self.get_mollist(m)
                    for n in mollist:
                        atomlist += self.mol.mols[n]
                elems = []
                for a in atomlist:
                    elems.append(self.mol.elems[a])
                natoms = len(elems)
            else:
                elems = self.get_elements()
                atomlist = range(self.natoms)
                natoms = self.natoms
            elems = [elem.upper() for elem in elems]
            xyz = self.get_xyz()
            xyz -= numpy.min(xyz, axis=0)
            fxyz = open(self.rundir + "/" + fname, "w")
            if zupfrac:
                zmin = min(xyz[:, 2])
                zmax = max(xyz[:, 2])
                delta_z = zmax - zmin
                linfo = []
                natoms = 0
                for i in xrange(len(mollist)):
                    minmol = min(self.mol.mols[i])
                    maxmol = max(self.mol.mols[i])
                    zsmol = xyz[minmol: maxmol + 1][:, 2]
                    if min(zsmol - zmin) >= 0 and (max(zsmol - zmin) <= zupfrac * delta_z):
                        linfo += [i]
                        natoms += len(self.mol.mols[i])
                natoms_less = 0
                for i in linfo:
                    for j in self.mol.mols[i]:
                        if elems[j] in rm_spc:
                            natoms_less += 1
                fxyz.write("%d\n" % (natoms - natoms_less))
                cell[2, 2] = zupfrac * delta_z
                fxyz.write(extxyz_head() + "\n")
                for i in linfo:
                    for j in self.mol.mols[i]:
                        if elems[j] not in rm_spc:
                            fxyz.write("%s %12.6f %12.6f %12.6f\n" % (elems[j], xyz[j, 0], xyz[j, 1], xyz[j, 2]))
            elif nmols:
                zmins = []
                for i in xrange(len(mollist)):
                    minmol = min(self.mol.mols[i])
                    maxmol = max(self.mol.mols[i])
                    zsmol = []
                    for j in range(minmol, maxmol + 1):
                        if elems[j] not in rm_spc:
                            zsmol += [xyz[j, 2]]
                    zmins += [numpy.average(zsmol)]
                indices = numpy.argsort(numpy.asarray(zmins))
                linfo = []
                natoms = 0
                for i in xrange(nmols):
                    linfo += [indices[i]]
                    natoms += len(self.mol.mols[indices[i]])
                natoms_less = 0
                for i in linfo:
                    for j in self.mol.mols[i]:
                        if elems[j] in rm_spc:
                            natoms_less += 1
                fxyz.write("%d\n" % (natoms - natoms_less))
                cell[2, 2] = max(xyz[self.mol.mols[indices[nmols]]][:, 2])
                fxyz.write(extxyz_head() + "\n")
                for i in linfo:
                    for j in self.mol.mols[i]:
                        if elems[j] not in rm_spc:
                            fxyz.write("%s %12.6f %12.6f %12.6f\n" % (elems[j], xyz[j, 0], xyz[j, 1], xyz[j, 2]))
            else:
                natoms_less = 0
                for i in xrange(natoms):
                    if elems[i] in rm_spc:
                        natoms_less += 1
                fxyz.write("%d\n" % (natoms - natoms_less))
                fxyz.write(extxyz_head() + "\n")
                for i in xrange(natoms):
                    if elems[i] not in rm_spc:
                        fxyz.write("%s %12.6f %12.6f %12.6f\n" % (elems[i], xyz[i, 0], xyz[i, 1], xyz[i, 2]))
            fxyz.close()
        return

    def write_tinker_xyz(self, fname, fullcell=None, mode="w", moldenr=False, img=False):
        """ if fullcell is not set we write the cell parameters (a,b,c,alpha,beta,gamma) and rotate the system properly
             if fullcell is true we write the complete cell vectors (9 values) to the first line

             :Parameters:
                 - moldenr : replaces atom type strings with integers to ensure readability by molden when set to True
             """
        self.pprint("   writing current configuration as tinker xyz file %s" % fname)
        cell = self.get_cell()
        if cell is not None:
            cellparams = unit_cell.abc_from_vectors(cell)
        if img:
            xyz = self.get_xyz_img()
        else:
            xyz  = self.get_xyz()
        if cell is not None:
            if fullcell:
                self.pprint("    writing complete cellvectors to file")
            else:
                rotcell = unit_cell.vectors_from_abc(cellparams)
                self.pprint("rotating system from current cellvectors\n")
                self.pprint(numpy.array2string(cell,precision=4,suppress_small=True))
                self.pprint("to new aligned cellvectors\n")
                self.pprint(numpy.array2string(rotcell,precision=4,suppress_small=True))
                # compute fractional coords in old cell
                inv_cell = numpy.linalg.inv(cell)
                abc = numpy.dot(xyz, inv_cell)
                # compute cartesian coordinates in new cellvectors from fractional coords
                xyz = numpy.dot(abc,rotcell)
        if moldenr:
            moltypes = {}
            for i, item in enumerate(self.mol.types):
                if not item in moltypes:
                    moltypes[item]=(len(moltypes)+1, self.mol.elems[i])
                self.mol.types[i] = moltypes[self.mol.types[i]][0]
        if self.is_master:
            f = open(self.rundir+"/"+fname, mode)
            if cell is not None:
                if fullcell:
                    f.write(("%5d # "+9*"%10.4f "+"\n") % tuple([self.natoms]+cell.ravel().tolist()))
                else:
                    f.write("%5d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n" % tuple([self.natoms]+cellparams))
            else:
                f.write("%5d \n" % self.natoms)
            for i in xrange(self.natoms):
                line = ("%3d %-3s" + 3*"%12.6f" + " %5s") % \
                    tuple([i+1]+[self.mol.elems[i]]+ xyz[i].tolist() + [self.mol.types[i]])
                try:
                    conn = (numpy.array(self.mol.cnct_orig[i])+1).tolist()
                except IndexError:
                    conn = (numpy.array(self.mol.cnct[i])+1).tolist()
                if len(conn) != 0:
                    line += (len(conn)*"%7d") % tuple(conn)
                f.write("%s \n" % line)
            if moldenr:
                f.write("#atomtypes \n")
                for i in moltypes:
                    f.write("%s %s %s \n" % (moltypes[i][0], i, moltypes[i][1]))
            f.close()
        return

    def write_config(self, fname, level=0):
        self.pprint ("  writing current configuration as a dl_poly CONFIG file to %s" % fname)
        self.dlp.write_config(fname, level)
        return


    def report_energies(self, full=False, marker=None, header_after_full=True):
        e = self.get_energy_contribs()
        if full:
            en = self.get_energy_names()
            self.pprint("Reporting current energies")
            eformat = "%%10s: %%%d.%df" % (self.eprec*2, self.eprec)
            for i,e in enumerate(e):
                self.pprint(eformat % (en[i], e))
            if header_after_full:
                # print the energy columns header after a full printout for readability
                self.report_energy_header(marker=marker)
        else:
            eformat = len(e)*("%%%d.%df " % (self.eprec*2+self.eprec_incr, self.eprec))
            line = eformat % tuple(e)
            if marker: line = marker+line
            self.pprint(line)
        return


    def report_energy_header(self, marker=None):
        en = self.get_energy_names()
        enformat = len(en)*("%%%ds " % (2*self.eprec+self.eprec_incr))
        line = enformat % tuple(en)
        if marker: line = marker+line
        self.pprint(line)
        return


    def get_energy_names(self):
        enames = ["vdW", "Coulomb", "bond", "angle", "oop", "torsion"]
        if self.use_efield: enames.append("efield")
        if self.QMMM: enames.append("QM")
        if self.extra_term: enames.append("extra")
        if self.plumed: enames.append("meta")
        for s in self.extra_systems:
            ex_enames = s.get_ener_names()
            if type(ex_enames) == types.ListType:
                enames += ex_enames
            else:
                enames.append(ex_enames)
        enames.append("TOTAL")
        if self.mode=="MD":
            enames.append("Ekin")
            for s in self.extra_systems:
                enames.append("Ek_"+s.get_name())
                enames.append("T_"+s.get_name())
            enames += ["Econs", "T_inst", "p_inst", "V", "X_T", "X_p"]
        return enames


    def get_energy_contribs(self):
        e = []
        e.append(self.dlp.engsrp/self.dlp.engunit)
        e.append(self.dlp.engcpe/self.dlp.engunit)
        e.append(self.dlp.engbnd/self.dlp.engunit)
        e.append(self.dlp.engang/self.dlp.engunit)
        e.append(self.dlp.enginv/self.dlp.engunit)
        e.append(self.dlp.engdih/self.dlp.engunit)
        if self.use_efield:
            e.append(self.dlp.engfld/self.dlp.engunit)
        etot = self.dlp.engcfg/self.dlp.engunit
        if self.QMMM:
            e.append(self.QMMM_eng)
            etot += self.QMMM_eng
        if self.extra_term:
            e.append(self.extra_eng)
            etot += self.extra_eng
        if self.plumed:
            e.append(self.meta_eng)
            etot += self.meta_eng
        for s in self.extra_systems:
            esys = s.get_epot()
            e.append(esys[0])
            e += esys[1]
            etot += esys[0]
        e.append(etot)
        econs = copy.copy(etot)
        if self.mode=="MD":
            ekin = self.dlp.engke/self.dlp.engunit
            ekin_rot = self.dlp.engrot/self.dlp.engunit
            ekin += ekin_rot
            e.append(ekin)
            econs += ekin
            for s in self.extra_systems:
                ekinsys = s.get_ekin()
                e.append(ekinsys)
                econs += ekinsys
                tempsys = self.calc_systemp(s)
                e.append(tempsys)
            #e.append(dlp.stpcns/self.dlp.engunit)
            econs += self.dlp.consv/self.dlp.engunit
            e.append(econs)
            e.append(self.dlp.stptmp)
            e.append(self.dlp.stpprs)
            e.append(self.dlp.stpvol)
            e.append(self.dlp.chit)
            e.append(self.dlp.chip)
        return numpy.array(e,"float64")

    # this set of methods is  a bit redundant to the above but in order not to break anything 
    # it is implemented like this to access scalar properties as a data_func
    scal_prop=[
            "vdW",     
            "Coulomb",
            "bond",  
            "angle",    
            "oop",   
            "torsion",  
            "efield",  
            "QM",
            "extra",
            "meta",
            "Ekin",
            "Ekin_rot",
            "T_inst",
            "p_inst",
            "V",
            "X_T",   
            "X_p",
        ]
        
    def get_scal_prop(self, propname):
        assert propname in self.scal_prop, "Unknown scalar property %s" % propname
        if propname=="vdW":
            x = self.dlp.engsrp/self.dlp.engunit
        elif propname=="Coulomb":
            x = self.dlp.engcpe/self.dlp.engunit
        elif propname=="bond":
            x = self.dlp.engbnd/self.dlp.engunit
        elif propname=="angle":
            x = self.dlp.engang/self.dlp.engunit
        elif propname=="oop":
            x = self.dlp.enginv/self.dlp.engunit
        elif propname=="torsion":
            x = self.dlp.engdih/self.dlp.engunit
        elif propname=="efield":
            if self.use_efield:
                x = self.dlp.engfld/self.dlp.engunit
            else:
                x = 0.0
        elif propname=="QM":
            if self.QMMM:
                x = self.QMMM_eng
            else:
                x = 0.0
        elif propname=="extra":
            if self.extra_term:
                x = self.extra_eng
            else:
                x = 0.0
        elif propname=="meta":
            if self.plumed:
                x = self.meta_eng
            else:
                x = 0.0
        elif propname=="Ekin":
            x = self.dlp.engke/self.dlp.engunit
        elif propname=="Ekin_rot":
            x = self.dlp.engrot/self.dlp.engunit
        elif propname=="T_inst":
            x = self.dlp.stptmp
        elif propname=="p_inst":
            x = self.dlp.stpprs
        elif propname=="V":
            x = self.dlp.stpvol
        elif propname=="X_T":
            x = self.dlp.chit
        elif propname=="X_p":
            x = self.dlp.chip
        else:
            # this should never happen
            raise ValueError
        return numpy.array(x)
        

    def report_timing(self):
        self.pprint ("\n\n*******************************************************")
        self.pprint ("           Timing report (all times in seconds)")
        # note that the names of all timers are hardcoded, so this might need to
        # be changed if you change the code (see force_module)
        timer_names=["startup","pair forces","three body", "four body", "bond/angle/torsion/inversion",\
                     "parallel force commun","qeq coulomb"]
        tn_length = 0
        # time_since_init = dlp_time.timer_since_init()
        for tn in timer_names:
            if len(tn)> tn_length: tn_length=len(tn)
        tn_format = "%%%ds   " % tn_length
        timers_percall = self.dlp_time.timers/self.dlp_time.tcounts
        total_time = numpy.sum(self.dlp_time.timers)
        self.pprint((tn_format%"timer")+("%12s  %12s  %8s"%("total","per call","ncalls")))
        for i in xrange(self.dlp_time.ntimers):
            if self.dlp_time.timers[i] != 0.0:
                line = tn_format  % timer_names[i]
                line +=  "%12.6f  "% self.dlp_time.timers[i]
                line +=  "%12.6f  "% self.timers_percall[i]
                line +=  "%8d"     % self.dlp_time.tcounts[i]
                self.pprint(line)
        self.pprint ((tn_length+3)*" "+("%12.6f total time"%total_time))
        # self.pprint ((tn_length+3)*" "+("%12.6f time since init"%time_since_init))
        self.pprint ((tn_length+3)*" "+("%12.6f cpu time"%time.clock()))
        # now we flush the timers and reset
        self.pprint ("Resetting Timers")
        self.dlp_time.timer_init()
        return

#################DEBUG things --- use only if you know what you are doing :-))

    def update_vdw_tables(self):
        self.dlp.update_vdw_tables()
        self.set_atoms_moved()
        return

    def write_vdw_pots(self, fname, write_grad=False, write_rep=False, add_disp=False):
        self.pprint("Writing vdw potential to file %s" % fname)
        if self.is_master:
            f = open(self.rundir+"/"+fname,"w")
            if write_grad:
                vdw = self.dlp_vdw.ggg
            else:
                vdw = self.dlp_vdw.vvv
            if write_rep:
                if add_disp:
                    vdw = self.dlp_vdw.vvv + self.dlp_vdw.vvvr
                else:
                    vdw = self.dlp_vdw.vvvr
            mxgrid = vdw.shape[0]
            mxvdw  = vdw.shape[1]-1
            for i in xrange(mxgrid):
                rrr = (i+1)*self.dlp.dlrpot
                line = "%10.5f" % rrr
                if write_grad:
                    line += (mxvdw*" %12.6f") % tuple(vdw[i,:-1]/(rrr*self.dlp.engunit))
                else:
                    line += (mxvdw*" %12.6f") % tuple(vdw[i,:-1]/self.dlp.engunit)
                f.write(line + "\n")
            f.close()
        return

    def report_nlist(self):
        self.pprint ("reporting complete neigborlist")
        self.pprint ("does not work in parallel runs!")
        for i in xrange(self.natoms):
            n = self.dlp_conf.lentry[i]
            line = (n*"%5d") % tuple(self.dlp_conf.list[i,:n])
            self.pprint(("%6d: "%(i+1))+line)
        return

    def report_nlistpar(self):
        self.pprint ("reporting complete neigborlist for parallel runs")
        ii = 0
        for i in xrange(self.idnode,self.natoms,self.nodes):
            n = self.dlp_conf.lentry[ii]
            line = (n*"%3d") % tuple(self.dlp_conf.list[ii,:n])
            # in order to see the data from all nodes we use print here ... output might be juggled
            print(("nlist node %1d %4d: "%(self.idnode,i+1))+line)
            ii += 1
        return

    def report_exlist(self):
        self.pprint ("reporting complete exclude list")
        self.pprint ("does not work in parallel runs!")
        for i in xrange(self.natoms):
            n = self.dlp_excl.nexatm[i]
            line = (n*"%4d") % tuple(self.dlp_excl.lexatm[i,:n])
            self.pprint(("%4d: "%(i+1))+line)
        return

    def report_exlistpar(self):
        self.pprint ("reporting complete exclude list for parallel runs")
        ii = 0
        for i in xrange(self.idnode,self.natoms,self.nodes):
            n = self.dlp_excl.nexatm[ii]
            line = (n*"%3d") % tuple(self.dlp_excl.lexatm[ii,:n])
            print(("exlist node %1d %4d: "%(self.idnode,i+1))+line)
            ii +=1
        return



############# MPI communication stuff - basically wrappers to basic_comm

    def bcast(self, data, type_is=None, dtype_is=None):
        if data != None:
            type_is = type(data)
            if type_is == numpy.ndarray:
                dtype_is = data.dtype
        if type_is == types.StringType:
            if self.is_master:
                n = len(data)
                self._pydlpoly.bcast_i(numpy.array([n],"i"),1)
                adata = numpy.array(data)
                self._pydlpoly.bcast_c(adata, n+1)
                # print self.is_master, data, len(data), type(data)
                return data
            else:
                n = numpy.array([0],"i")
                self._pydlpoly.bcast_i(n,1)
                n=n[0]
                adata=numpy.array((n+1)*" ")
                self._pydlpoly.bcast_c(adata, n+1)
                # please do not ask me why we need to chop off the last character ... i do not get it
                data= adata.tostring()[:-1]
                # print self.is_master, data, len(data), type(data)
                return data
        elif type_is == numpy.ndarray:
            # this is a numpy array
            if dtype_is == numpy.dtype("int32"):
                # its a 4 byte integer
                if self.is_master:
                    ndim = len(data.shape)
                    self._pydlpoly.bcast_i(numpy.array([ndim],"i"),1)
                    nshape = numpy.array(data.shape, "i")
                    self._pydlpoly.bcast_i(nshape,ndim)
                    total_len = numpy.product(nshape)
                    self._pydlpoly.bcast_i(data.ravel(),total_len)
                    return data
                else:
                    ndim = numpy.array([0],"i")
                    self._pydlpoly.bcast_i(ndim, 1)
                    ndim = ndim[0]
                    nshape = numpy.array(ndim*[0],"i")
                    self._pydlpoly.bcast_i(nshape, ndim)
                    data = numpy.zeros(nshape, "int32")
                    total_len = numpy.product(nshape)
                    self._pydlpoly.bcast_i(data.ravel(),total_len)
                    return data
            elif dtype_is == numpy.dtype("float64"):
                # its a 8 byte double float
                if self.is_master:
                    ndim = len(data.shape)
                    self._pydlpoly.bcast_i(numpy.array([ndim],"i"),1)
                    nshape = numpy.array(data.shape, "i")
                    self._pydlpoly.bcast_i(nshape,ndim)
                    total_len = numpy.product(nshape)
                    self._pydlpoly.bcast_d(data.ravel(),total_len)
                    return data
                else:
                    ndim = numpy.array([0],"i")
                    self._pydlpoly.bcast_i(ndim, 1)
                    ndim = ndim[0]
                    nshape = numpy.array(ndim*[0],"i")
                    self._pydlpoly.bcast_i(nshape, ndim)
                    data = numpy.zeros(nshape, "float64")
                    total_len = numpy.product(nshape)
                    self._pydlpoly.bcast_d(data.ravel(),total_len)
                    return data
            else:
                print "this numpy data type is not supported for bcast"
                return
        else:
            print "needs to be implemented"
            return

    def gdsum(self, data):
        #print "DEBUG USING GDSUM"
        #n = len(data)
        buf = numpy.zeros(data.shape, data.dtype)
        #self._pydlpoly.gdsum(data, n, buf)
        # use mpi4py instead
        self.local_comm.Allreduce(data,buf,MPI.SUM)
        data[::] = buf[::]
        return data

    ############# QEq temp #####################

    def init_qeq(self,use_Jij=False, input_file=None, exclude=None, acks2=False):
        self.dlp.qeq = True
        self.qeq = qeq.qeq(self, input_file=input_file, exclude=exclude, acks2=acks2)
        # this changes the energy expression ... enforce a recalc
        self.set_atoms_moved()
        # start up Jij storing ... maybe we should do this only upon request
        if use_Jij:
            self.dlp_coul.init_Jij(self.idnode)
            self.use_Jij=True
            self.valid_Jij=False
        return

    def delete_qeq(self):
        self.dlp.qeq = False
        self.dlp_setup.qeq_iter = False
        qeq = self.qeq
        self.qeq = None
        del(qeq)
        return



    ############# Extra System Methods ##########

    def calc_systemp(self, system):
        sysdof = system.get_dof()
        sysekin = system.get_ekin()
        systemp = 2*sysekin/(sysdof*boltzkcal)
        return systemp


    ####################### LBFGS ###############

