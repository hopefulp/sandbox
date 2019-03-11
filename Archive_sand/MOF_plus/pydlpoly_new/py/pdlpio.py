"""

                            pdlpio

    implements a hdf5/h5py based class to handle storage of pydlpoly
    molecular systems (used for restart and trajectory data)

    the class object is a frontend to a hdf5/h5py file, and it can be used
    outside of pydlpoly to read and postprocess trajectory data

    one part of the data is "fixed"
    - number of atoms
    - elements
    - atom types
    - connectivity
    - boundary conditions
    this part is defined from the tinker xyz input file
    in a restart this data is used as a replacement to the usual tinker xyz input

    the other part is "staged"
    in other words if needed multiple stages can be defined like "equilibration",
    "sample1", "sample2" etc.
    the initial and default stage is called "default" by default ... yes!
    for each stage current restart information (position/velocitites/cell) is
    automatically stored every N steps. restart data is stored in double precision
    in addition, trajectory info is stored every M steps. (both N and M can be changed)
    trajectory data (pos or pos/vel or pos/vel/force) can be stored both in single and
    double precision on request.
    
    NOTE: This is NOT a paralle code => in pydlpoly we need to make sure that
          all routines are called by the master only
          
          We should chnage this in the future. currently on restart a fiel si opened on ALL nodes
          and the reading is done parallel (as in ascii reading from xyz in assign_FF).
          better: pass a comm (local_comm of pydlpoly) to pdlpio and use it for exchange.
          read by default also only on the master and broadcast objects to slaves.
          in case of write ony write on master. from the caller all is done in parallel.
          but within the object only the master has an open hdf5 file
    
    Versions (in file attribute "version"
    - no version info == 1.0 : original
    - version 1.1            : nstep of trajectories is stored per dataset and not for the stage
    
"""

import h5py
import numpy

import types

class pdlpio:

    def __init__(self, fname, mode="a"):
        """ open the hdf5 file depending on the mode
            defaults are set
            
            NEW 2014: if the opened file exists 8and contains data read it in
                      this is used in the pdlp script (not with pydlpoly)
        """
        self.verbose = 0
        #
        self.fname = fname
        self.mode  = mode
        #
        self.h5file = h5py.File(fname, mode)
        #
        self.file_version = 1.1
        
        #
        if "version" in self.h5file.attrs.keys():
            if (self.mode == "a") or (self.mode == "w"): 
                assert self.file_version == self.h5file.attrs["version"], "Exisiting file has a different version! Can not add data"
        else:
            self.h5file.attrs["version"] = self.file_version
        # defaults
        self.pd = None
        self.stagelist     = []
        self.track_data = None
        self.traj_nstep = 1
        self.rest_nstep = 1
        #
        # 
        if "system" in self.h5file.keys():
            # ok there is some system so initalize data from here
            self.stagelist = self.h5file.keys()
            self.stagelist.remove("system")
            self.system = self.h5file["system"]
            self.natoms = self.system["elems"].shape[0]
            self.bcd    = self.system.attrs["bcd"]
        else:
            self.system    = self.h5file.require_group("system")
            self.natoms = 0
            self.bcd    = 0
        #
        # helper object for hdf5 variable length strings
        self.str_dt = h5py.new_vlen(str)
        # track charges if floating_charges is True
        self.floating_charges = False
        return
        
    def pprint(self, it):
        print(it)
        return

    def close(self):
        self.h5file.close()
        return
        
    def initialise(self, initial_stage="default", restart_steps=10, pd=None):
        # if self.pd == None then this is a run without attached pydlpoly -> analysis (can' write restart)
        self.pd = pd
        if self.pd:
            self.data_funcs = self.pd.get_data_funcs()
            # go over all data objects and register their shapes
            self.data_shapes = {}
            for d in self.data_funcs.keys():
                self.data_shapes[d] = self.data_funcs[d]().shape
            self.rest_data = ["xyz", "vel"]
            if self.bcd > 0:
                self.rest_data.append("cell")
            if self.floating_charges:
                self.rest_data.append("charges")
            # for the IO get a ref to pydlpolys pprint
            self.pprint = self.pd.pprint
        self.add_stage(initial_stage, restart_steps)
        self.counter = 0
        return
        
    # system information
        
    def set_system(self, elems, atypes, cnc_table, bcd):
        """ NOTE: cnc_table must be a shape = (nbonds, 2) object """
        if self.natoms > 0:
            # NOTE: Add more checks here!
            consistent = True
            if len(elems) != self.natoms : consistent = False
            if not consistent: 
                raise IOError, "pdlp file already contains atom information"
        self.natoms = len(elems)
        self.empty = False
        self.elems = self.system.require_dataset("elems", shape=(self.natoms,), dtype=self.str_dt)
        self.elems[...] = elems
        self.atypes = self.system.require_dataset("atypes", shape=(self.natoms,), dtype=self.str_dt)
        self.atypes[...] = atypes
        if len(cnc_table) > 0:
            # catch if there are no bonds at all: h5py does not like zero size selections
            cnc_table = numpy.array(cnc_table, dtype="i")
            self.cnc_table = self.system.require_dataset("cnc_table", shape=cnc_table.shape, dtype=cnc_table.dtype )
            self.cnc_table[...] = cnc_table
        self.system.attrs["bcd"] = bcd
        self.bcd = bcd
        return

    def get_system(self):
        # NOTE: check that copys of h5py objects are passed (should be ok with tolist())
        elems     =  list(self.system["elems"])
        atypes    =  list(self.system["atypes"])
        if "cnc_table" in self.system.keys():
            cnc_table = numpy.array(self.system["cnc_table"]).tolist()
        else: 
            cnc_table = []
        bcd       =  self.system.attrs["bcd"]
        return elems, atypes, bcd, cnc_table
        
    def set_molecules(self, whichmol, moltypes, molnames):
        """ whichmol:  integer array of length natoms
            moltypes:  integer array of length nmols
            molnames:  list of strings of the length of nmoltypes
        """
        if not self.system.keys().count("elems"):
            raise IOError, "pdlp file does not contain any system info, can't write molecuels data"
        nmols     = len(moltypes)
        nmoltypes = len(molnames)
        self.whichmol = self.system.require_dataset("whichmol", shape=(self.natoms,), dtype="i")
        self.whichmol[...] = whichmol
        if "moltypes" in self.system.keys():
            # entries exist already ... resize if necessary
            if self.system["moltypes"].shape[0] != nmols:
                self.system["moltypes"].resize((nmols,))
            if self.system["molnames"].shape[0] != nmoltypes:
                self.system["molnames"].resize((nmoltypes,))
        self.moltypes = self.system.require_dataset("moltypes", shape=(nmols,), maxshape=(None,), dtype="i")
        self.moltypes[...] = moltypes
        self.molnames = self.system.require_dataset("molnames", shape=(nmoltypes,), maxshape=(None,), dtype=self.str_dt)
        self.molnames[...] = molnames
        return

    def get_molecules(self):
        whichmol = list(self.system["whichmol"])
        moltypes = list(self.system["moltypes"])
        molnames = list(self.system["molnames"])
        return whichmol, moltypes, molnames

    def has_stage(self, stagename):
        if self.h5file.keys().count(stagename) > 0:
            return True
        else:
            return False
        
    def get_stages(self):
        return self.stagelist
        
    def add_stage(self, stagename, nstep):
        """ adds a new stage and initalises the groups and datasets needed """
        if stagename == "system":
            raise IOError, "A stage should never be named system!"
        self.rest_nstep= nstep
        self.stagename = stagename
        self.stagelist.append(stagename)
        # add the restart part
        self.stage     = self.h5file.require_group(self.stagename)
        self.restart   = self.stage.require_group("restart")
        if self.pd:
            self.rest_datasets = {}
            for d in self.rest_data:
                self.rest_datasets[d] = self.restart.require_dataset(d, shape=self.data_shapes[d], dtype="float64")
        self.counter = 0
        return
        
    def add_traj(self, nstep, data, prec="float32", tstep=0.001):
        """ adds a trajectory group to the current stage and starts to track data (list of names)
            refering to the directories data_funcs and data_shapes 
            
            :Parameters:
                - nstep (int or list of int) : number of steps before a traj data is written. If it is a list the info is per
                                 entry in the data list, otherwise for all datasets
                - data (list of strings) : string keys which should be written to the pdlp file
                - prec (numpy type specifier) : default = "float32"
                - tstep (float) : timestep in ps
        """
        if self.stage.keys().count("traj"):
            self.pprint("This stage already tracks data! Please generate a new stage first")
            return
        # add a group for trajectory information
        self.traj      = self.stage.require_group("traj")
        self.traj.attrs["nstep"] = 0
        self.traj.attrs["tstep"] = tstep
        if type(nstep) != types.ListType:               #RS NOTE: this is a list now
            self.traj_nstep = len(data)* [nstep]
        else:
            self.traj_nstep = nstep
        self.traj_frame = len(self.traj_nstep) * [0]
        self.track_data = data
        # generate data sets with extendable first dimension
        self.traj_datasets = {}
        for i,d in enumerate(self.track_data):
            if not self.data_funcs.has_key(d):
                raise IOError, "The data obejct %s is unknown!" % d
            dshape = list(self.data_shapes[d])
            self.traj_datasets[d] = self.traj.create_dataset(d,dtype=prec,\
                                       shape=tuple([0]+dshape),\
                                       maxshape=tuple([None]+dshape),\
                                       chunks=tuple([1]+dshape))
            self.traj_datasets[d].attrs["nstep"] = self.traj_nstep[i]
        return
        
    def __call__(self, force_wrest=False):
        """ this generic routine is called to save restart info to the hdf5 file
            if force_wrest is True then restart is written in any case """
        if not self.pd: raise IOError, "No pydlpoly instance"
        data_written = False
        if (self.counter%self.rest_nstep == 0) or force_wrest:
            # restart is written .. fetch all data 
            for d in self.rest_data:
                if self.verbose:
                    self.pprint("Writing restart data %s to pdlp file" % d)
                data = self.data_funcs[d]()
                self.rest_datasets[d][...] = data
            data_written = True
        if self.track_data != None:
            for i,d in enumerate(self.track_data):
                if (self.counter%self.traj_nstep[i] == 0):
                    self.traj_frame[i] += 1
                    data = self.data_funcs[d]()
                    tds = self.traj_datasets[d]
                    tds.resize(self.traj_frame[i], axis=0)
                    tds[self.traj_frame[i]-1,...] = data
                    data_written = True
        # now we are done 
        if data_written : self.h5file.flush()
        self.counter += 1
        return
        
    def read_restart(self, stage, data):
        """ reads data from restart stage and returns it """
        stage = self.h5file[stage]
        restart = stage["restart"]
        data = restart[data]
        return numpy.array(data)
            
        