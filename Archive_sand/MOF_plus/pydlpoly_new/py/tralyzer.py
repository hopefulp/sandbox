""" 
               tralyzer
               
       A trajectory analyzer working on and processing pdlp files 
       
"""

from mpi4py import MPI
import numpy as nm
import numpy.linalg as la
import pdlpio
import copy

# NOTE currently pydlpoly uses masses from the key file!! these are strange
#      correct them or better use instead a common elements.py (also in weaver and assign_FF)
atomicmass = {\
    "h" : 1.0079 , "he": 4.0026,\
    "b" : 10.811, "c" : 12.011, "n": 14.007, "o": 15.999, "f": 19.0,\
    "ne": 20.0, "s": 32, "cl": 35, "ar":40.0, "si": 28.0855, "ga":69.723, "zn": 65.409, \
    "br": 79.9,"kr": 83.0, "xe":131.0 , "cu":63.546, "se": 90.0, "zr": 91.22\
    }
atomicnumbers = {"xx" : 1,\
    "h" : 1, "he" : 2, \
    "c" : 6, "n"  : 7, "o" : 8,  "al": 13, \
    "si": 14, "cu": 29, "zn": 30, "ga":31, "ge": 32, "hf": 72}
a2b = 1.0/0.529177

class tralyzer:
    
    def __init__(self, fname):
        """ start up and read what can be read from pdpio file 
        NOTE: We should add some nicer printout  """
        self.comm = MPI.COMM_WORLD
        self.node_id = self.comm.Get_rank()
        self.nodes   = self.comm.Get_size()
        # open file
        self.pdpf = pdlpio.pdlpio(fname, mode="a")
        self.version = 1.0
        if "version" in self.pdpf.h5file.attrs.keys():
            self.version = self.pdpf.h5file.attrs["version"]
        self.elems, self.types, self.boundarycond, ctab = self.pdpf.get_system()
        self.natoms = len(self.elems)
        # recover connectivity from ctab
        self.cnct = []
        for i in xrange(self.natoms): self.cnct.append([])
        for c in ctab:
            i,j = c
            self.cnct[i].append(j)
            self.cnct[j].append(i)
        # register types
        self.typedata = {}
        for t in self.types:
            if not self.typedata.has_key(t): self.typedata[t] = None
        self.ntypes = len(self.typedata.keys())
        # molecule info - note that only a nonredundant set of info is in pdlp and the
        # rest needs to be recovered
        self.whichmol, self.moltypes, self.molnames = self.pdpf.get_molecules()
        self.nmols = len(self.moltypes)
        self.mols = []
        for i in xrange(self.nmols) : self.mols.append([])
        for i, m in enumerate(self.whichmol): self.mols[m].append(i)
        # REPORT
        self.pprint(80*"#")
        self.pprint("\n")
        self.pprint(30*" "+"TRALYZER")
        self.pprint("\n")
        self.pprint(80*"#")
        self.pprint("opened pdlp file with %d atoms and %d molecules" % (self.natoms, self.nmols))
        if self.nodes>1:
            self.pprint("running parallel on %d nodes, distributing work" % self.nodes)
        return
                 
    def pprint(self, string):
        if self.node_id==0:
            print(string)
        return
        
    def set_active_stage(self, stagename, timestep=0.001, use_fixed_cell=False):
        if not self.pdpf.has_stage(stagename):
            raise IOError, "The pdp file has no stage %s" % stagename
        self.stagename = stagename
        self.traj = self.pdpf.h5file[stagename]["traj"]
        self.rest = self.pdpf.h5file[stagename]["restart"]
        # for version > 1.0 allow different steps
        if self.version > 1.0:
            nstep = self.traj.attrs["nstep"]
            if len(nstep) > 1:
                nstep = nstep.max()
            if "nstep" in self.traj["xyz"].attrs.keys():
                nstep = self.traj["xyz"].attrs["nstep"]
        self.dt = timestep*nstep
        #self.xyz = self.traj["xyz"]
        # self.nframes = self.xyz.shape[0]
        self.nframes = self.traj["xyz"].shape[0]
        self.pprint("set active stage to %s" % stagename)
        self.pprint("number of frames    %10d" % self.nframes)
        self.pprint("timestep            %10.5f" % self.dt)
        if self.boundarycond > 0:
            if (not self.traj.keys().count("cell")) or use_fixed_cell:
                self.pprint("using static cellparams from restart data")
                self.cell = nm.array(self.rest["cell"])
                self.fixed_cell = True
            else:
                self.pprint("using trajectory cellparams per frame")
                self.fixed_cell = False
                self.cell = self.traj["cell"]
        if self.traj.keys().count("vel"):
            self.pprint("velocities are available")
            self.vel = self.traj["vel"]
        else:
            self.vel = None
        return
        
    def get_mollist(self, molname):
        moltype = self.molnames.index(molname)
        return [i for i,n in enumerate(self.moltypes) if n==moltype]
        
    def write_xyz(self, fname, mols=None, first_frame=None, last_frame=None, stride=None, add_com=None):
        """ writes a simple xyz file out
             mols        : list of molnames if you want to restrict to these (default None)
             first_frame : first frame to write (default = 0)
             last_frame  : last frame to write  (default = nframes)
             stride      : default = 1
                           REMARK we sue python index nomenclature!
             add_com     : a tuple of the type ("molname", "element")
                           adds an atom of the type "element" to the xyz file at the COm position stored in 
                           dataset molname_COM (must be stored!!)
        """
        if self.node_id== 0:
            if mols:
                atomlist = []
                for m in mols:
                    mollist = self.get_mollist(m)
                    for n in mollist: atomlist += self.mols[n]
                elems = []
                for a in atomlist: elems.append(self.elems[a])
                natoms = len(elems)
            else:
                elems = self.elems
                atomlist = range(self.natoms)
                natoms = self.natoms
            if add_com:
                com = self.load_com(add_com[0])
                com_elem = add_com[1]
                ncom = com.shape[1]
                natoms += ncom
            fxyz = open(fname, "w")
            if not first_frame: first_frame = 0
            if not last_frame:  last_frame  = self.nframes
            if not stride:      stride = 1
            for f in xrange(first_frame, last_frame, stride):
                fxyz.write("%10d\n\n" % natoms)
                cxyz = self.xyz[f]
                for i,a in enumerate(atomlist):
                    pos = cxyz[a]
                    fxyz.write("%4s %12.6f %12.6f %12.6f\n" % (elems[i], pos[0], pos[1], pos[2]))
                if add_com:
                    ccom = com[f]
                    for i in xrange(ncom):
                        pos = ccom[i]
                        fxyz.write("%4s %12.6f %12.6f %12.6f\n" % (com_elem, pos[0], pos[1], pos[2]))
            # done
        return
             
    ####  cell stuff for trajectories of NsT runs
        
    def get_cell_averages(self, first_frame=0, last_frame=None, write_frames=None):
        """ computes averages over cell params and volume """
        if self.fixed_cell:
            self.pprint("This makes sense only for trajectories with varying cell")
            return
        cell = nm.array(self.cell[first_frame:last_frame], dtype="float64")
        V       = nm.sum(cell[:,0]*nm.cross(cell[:,1],cell[:,2]), axis=1)
        mean_V  = nm.mean(V)
        stdd_V  = nm.std(V)
        cellpar = nm.sqrt(nm.sum(cell*cell, axis=2))
        mean_cp = nm.mean(cellpar, axis=0)
        stdd_cp = nm.std(cellpar, axis=0)
        alpha   = nm.arccos(nm.sum(cell[:,1]*cell[:,2], axis=1)/(cellpar[:,1]*cellpar[:,2]))*(180.0/nm.pi)
        beta    = nm.arccos(nm.sum(cell[:,0]*cell[:,2], axis=1)/(cellpar[:,0]*cellpar[:,2]))*(180.0/nm.pi)
        gamma   = nm.arccos(nm.sum(cell[:,0]*cell[:,1], axis=1)/(cellpar[:,0]*cellpar[:,1]))*(180.0/nm.pi)
        mean_alp= nm.mean(alpha)
        stdd_alp= nm.std(alpha)
        mean_bet= nm.mean(beta)
        stdd_bet= nm.std(beta)
        mean_gam= nm.mean(gamma)
        stdd_gam= nm.std(gamma)
        self.pprint("Analysis of cellparameters")
        self.pprint("==========================")
        self.pprint("Volume:   %12.6f A**3  (std-dev %12.6f)" % (mean_V, stdd_V))
        self.pprint("a     :   %12.6f A     (std-dev %12.6f)" % (mean_cp[0], stdd_cp[0]))
        self.pprint("b     :   %12.6f A     (std-dev %12.6f)" % (mean_cp[1], stdd_cp[1]))
        self.pprint("c     :   %12.6f A     (std-dev %12.6f)" % (mean_cp[2], stdd_cp[2]))
        self.pprint("alpha :   %12.6f deg   (std-dev %12.6f)" % (mean_alp, stdd_alp))
        self.pprint("beta  :   %12.6f deg   (std-dev %12.6f)" % (mean_bet, stdd_bet))
        self.pprint("gamma :   %12.6f deg   (std-dev %12.6f)" % (mean_gam, stdd_gam))
        if write_frames:
            self.pprint("Writing results for each frame to file %s" % write_frames)
            of = open(write_frames, "w")
            for i in xrange(len(cell)):
                of.write("%5d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n" %\
                         (i+1, cellpar[i,0], cellpar[i,1], cellpar[i,2], alpha[i], beta[i], gamma[i], V[i]))
            of.close()
        return
        
    def write_cell_as_list(self, fname):
        h = nm.array(self.cell)
        N = h.shape[0]
        h = nm.reshape(h, (N, 9))
        f = open(fname, "w")
        for i in xrange(N):
            f.write(("%5d" + 9*"%12.6f " + "\n") % tuple([i]+h[i].tolist()))
        f.close()
        return
        
        
    def compute_elastic_constants(self, T, first_frame=0, last_frame=-1, h0=None):
        """ Compute elastic constants by strain fluctuation .. code basis from FX Coudert """
        if self.fixed_cell:
            self.pprint("This makes sense only for trajectories with varying cell")
            return
        h = nm.array(self.cell[first_frame:last_frame])
        # reference cell h0 (inverted and transposed) ... could be the equl cell or anything else
        if not h0: h0 = h[0]
        h0m1  = la.inv(h0)
        h0m1t = h0m1.transpose()
        # compute strain (deviation from reference cell)
        eps_mat = nm.array(map(lambda hi : (nm.dot(h0m1t, nm.dot(hi.transpose(), nm.dot(hi, h0m1))) - nm.identity(3))/2.0, h))
        # Voigt_map = [[0,0], [1,1], [2,2], [2, 1], [2, 0], [1, 0]]
        eps_voigt = nm.empty([6,eps_mat.shape[0]], "d")
        # eps_voigt_std = nm.empty([6,eps_mat.shape[0]], "d")
        eps_voigt[0] = eps_mat[:,0,0]
        eps_voigt[1] = eps_mat[:,1,1]
        eps_voigt[2] = eps_mat[:,2,2]
        eps_voigt[3] = eps_mat[:,2,1]
        eps_voigt[4] = eps_mat[:,2,0]
        eps_voigt[5] = eps_mat[:,1,0]
        # DEBUG
        # print eps_mat[0:10]
        
        # DEBUG
        eps_voigt_mean = nm.mean(eps_voigt, axis=1)
        # eps_voigt_std = nm.std(eps_voigt, axis=1)
        Smat = nm.zeros([6,6], "d")
        # Smat_std = nm.zeros([6,6], "d")
        for i in xrange(6):
            for j in xrange(i,6):
                Smat[i,j] = nm.mean(eps_voigt[i]*eps_voigt[j])-eps_voigt_mean[i]*eps_voigt_mean[j]
                if i != j: Smat[j,i] = Smat[i,j]
                # fix this for error propagation!!
                # Smat_std[i,j] = nm.std(eps_voigt[i]*eps_voigt[j], axis=1)-eps_voigt_std[i]*eps_voigt_std[j]
                # if i != j: Smat_std[j,i] = Smat_std[i,j]
        # finish it (use mean Volume and T)
        V = nm.mean(nm.sum(h[:,0]*nm.cross(h[:,1],h[:,2]), axis=1))
        Smat *= (V*1.e-30)/(1.3806488e-23*T)
        print "Smat raw"
        print Smat
        # Cmat in GPa
        Cmat = la.inv(Smat)*1.0e-9
        print "Cmat raw"
        print Cmat
        Cmat_eig = la.eigvals(Cmat)
        print Cmat_eig
        return
        
    ####  COM stuff
        
    def compute_com(self, molname, store=True):
        mollist = self.get_mollist(molname)
        nmols = len(mollist)
        self.pprint("\ncomputing the COM of %d molecules with name %s" % (nmols,molname))
        self.pprint("Molecules are unwrapped .. this works only if the molecule is smaller then half the unit cell!!")
        atomlist = []
        for m in mollist:
            atomlist.append(self.mols[m])
        # since all molecules are equal we need to get masses for only one representative
        mass = []
        for a in atomlist[0]:
            mass.append(atomicmass[self.elems[a]])
        # make it the right shape to multiply with xyz coords
        mass = nm.array(mass)[:,nm.newaxis]
        summass = nm.sum(mass)
        # allocate the array for the COMs
        # NOTE we could use the dtype of the xyz dataype here 
        com = nm.zeros([self.nframes,nmols,3],dtype=self.xyz.dtype)
        if self.fixed_cell:
            celldiag = self.cell.diagonal()
            half_celldiag = celldiag/2.0
        #if self.nodes>1:
            #frames_per_node = self.nframes/self.nodes
            #wrapup          = self.nframes%self.nodes
            #first_frame = self.node_id*frames_per_node
            #last_frame  = first_frame+frames_per_node
            #if self.node_id == self.nodes-1:
                #last_frame += wrapup
        #else:
            #first_frame=0
            #last_frame =self.nframes
        for f in xrange(self.node_id, self.nframes, self.nodes):
            if not self.fixed_cell:
                celldiag = self.cell[f].diagonal()
                half_celldiag = celldiag/2.0
            cxyz = self.xyz[f]
            for i,m in enumerate(atomlist):
                mxyz = cxyz[m]
                # unwrap the molecule, the reference image is the first atom:
                #         if any other atom is farther then the half cell we need to unwrap it
                #         NOTE this works only for molecules smaller then half the cell!!!
                dmxyz = mxyz[1:]-mxyz[0]
                mxyz[1:] += nm.where(nm.less_equal(dmxyz, -half_celldiag), 1.0, 0.0)*celldiag
                mxyz[1:] -= nm.where(nm.greater   (dmxyz,  half_celldiag), 1.0, 0.0)*celldiag
                ccom = nm.sum(mxyz*mass, axis=0)/summass
                # now wrap ccom back into the box if it is outside (box goes from -halfcell to +halfcell)
                ccom += nm.where(nm.less   (ccom,-half_celldiag), celldiag, 0.0)
                ccom -= nm.where(nm.greater(ccom, half_celldiag), celldiag, 0.0)
                com[f,i,:] = ccom
        if self.nodes>1:
            buf = nm.zeros([self.nframes,nmols,3],dtype=self.xyz.dtype)
            self.comm.Allreduce(com, buf, MPI.SUM)
            com = buf
        if store:
            dataname = molname+"_COM"
            self.pprint("storing COM info in dataset %s" % dataname)
            if self.node_id ==0:
                self.com_h5 = self.traj.require_dataset(dataname, shape=com.shape, dtype=com.dtype)
                self.com_h5[...] = com
        return com
        
    def load_com(self, molname):
        self.com_h5 = self.traj[molname+"_COM"]
        return nm.array(self.com_h5)
        
    #### 3D Maps (COM probability distributions)
        
    def map_com(self, com, gridsize, name, store=True, replace=False):
        self.pprint("\nComputing a 3D map of a com for molecule %s" % name)
        nframes = com.shape[0]
        gridsize = nm.array(gridsize,dtype="int32")
        grid = nm.zeros(gridsize+nm.array([1,1,1]), dtype="d")
        #grid = nm.zeros(gridsize, dtype="d")
        self.pprint("WARING: Map currently works only for constant volume ensembles!!")
        if not self.fixed_cell:
            raise ValueError, "map_com works currently only for fixed cell sizes"
        if not ((self.boundarycond==1)or(self.boundarycond==2)):
            raise ValueError, "Can't map for these boundaryconditions"
        # shift com by half a cell +h/2 up to get all coordinates above zeros
        celldiag = self.cell.diagonal()
        h = celldiag/gridsize
        com += (celldiag+h)/2.0
        int_com = (com/h).astype("int32")
        entries = 0
        for f in xrange(self.node_id, nframes, self.nodes):
            for ind in int_com[f].tolist():
                grid[tuple(ind)]+=1.0
                entries += 1
        # comunicate 
        if self.nodes > 1:
            # first allreduce (also number of entries need to be allreduced)
            buf = nm.zeros(grid.shape, grid.dtype)
            self.comm.Allreduce(grid, buf, MPI.SUM)
            grid = buf
            entries = self.comm.allreduce(entries)
        # fix boundary conditions (values at the "faces" need to be identical)
        grid[0,:,:] += grid[-1,:,:]
        grid[:,0,:] += grid[:,-1,:]
        grid[:,:,0] += grid[:,:,-1]
        grid[-1,:,:] = grid[0,:,:]
        grid[:,-1,:] = grid[:,0,:]
        grid[:,:,-1] = grid[:,:,0]
        # normalize
        grid /= (float(entries)*h.prod())
        if store:
            dataname = name+"_MAP"
            if self.nodes>1:
                if replace:
                    self.pprint("trying to replace the dataset %s" % dataname)
                    del(self.traj[dataname])
                self.pprint("storing map in dataset %s" % dataname)
                self.map_h5 = self.traj.require_dataset(dataname, shape=grid.shape, dtype=grid.dtype)
                self.map_h5[...] = grid
        return grid

    def load_map(self, mapname):
        self.map_h5 = self.traj[mapname]
        return nm.array(self.map_h5)
        
    def write_map_as_cube(self, map3d, fname, geomstage, showmol):
        if not ((self.boundarycond==1)or(self.boundarycond==2)):
            raise ValueError, "Can't map for these boundaryconditions"
        self.pprint("\nWriting cube file of map with geometry from stage %s" % geomstage)
        mollist = self.get_mollist(showmol)
        atomlist = []
        for m in mollist: atomlist += self.mols[m]
        elems = []
        for a in atomlist: elems.append(atomicnumbers[self.elems[a]])
        all_xyz = nm.array(self.pdpf.h5file[geomstage]["restart"]["xyz"])
        mol_xyz = all_xyz.take(atomlist, axis=0)
        natoms = len(atomlist)
        N3 = list(map3d.shape)
        if len(N3) == 2:
            self.pprint("Writing 2D map as a pseudo 3D")
            N3 += [1]
        N3 = nm.array(N3)
        celldiag = self.cell.diagonal()
        hcell = celldiag/2.0*a2b
        mol_xyz *= a2b
        h3 = (celldiag/N3)*a2b
        if self.node_id==0:
            f = open(fname+".cube","w")
            f.write("3DMAP as cube %s \n" % (fname))
            f.write("\n")
            f.write("%5d%12.6f%12.6f%12.6f\n" % (natoms, -hcell[0]+h3[0]/2.0, -hcell[1]+h3[1]/2.0, -hcell[2]+h3[2]/2.0))
            f.write("%5d%12.6f%12.6f%12.6f\n" % (N3[0], h3[0], 0.0, 0.0))
            f.write("%5d%12.6f%12.6f%12.6f\n" % (N3[1], 0.0, h3[1], 0.0))
            f.write("%5d%12.6f%12.6f%12.6f\n" % (N3[2], 0.0, 0.0, h3[2]))
            for i in xrange(natoms):
                xyz = mol_xyz[i]
                f.write("%5d%12.6f%12.6f%12.6f%12.6f\n" % (elems[i], 0.0, xyz[0], xyz[1], xyz[2]))
            lformat = N3[2]*"%e "+"\n"
            for x in xrange(N3[0]):
                for y in xrange(N3[1]):
                    ld = map3d[x,y]
                    if N3[2] == 1: ld = [ld]
                    f.write(lformat % tuple(ld))
            f.close()
        return

# NOT working properly ... why??
#    def map_symmetrize(self, map3d):
#        """ currently only a cubic symmetrization is performed """
#        new_map = (map3d+nm.transpose(map3d, axes=(1,2,0))+nm.transpose(map3d, axes=(2,0,1)))/3.0
#        return new_map

    def com_msd(self, com, molname, dt_max=None, store=True):
        """ computes the mean square diffusion for a given COM 
            we keep the MSD per molecule and per dimension
            averaging is easy by taking sum of the array """
        if not ((self.boundarycond==1)or(self.boundarycond==2)):
            raise ValueError, "Can't map for these boundaryconditions"
        self.pprint("\nComputing MSD from COM of molecule %s" % molname)
        # set up
        if dt_max:
            self.pprint("maximum deltat=%10.5f ps of a total of %10.5f ps" % (dt_max, self.nframes*self.dt))
            max_msd_frames = int(dt_max/self.dt)
            if max_msd_frames > self.nframes-1:
                raise ValueError, "requested detla t is more then the avialbale sampling time"
        else:
            max_msd_frames = self.nframes
        celldiag = self.cell.diagonal()
        celldiag_half = celldiag/2.0
        # before we can compute the MSD we have to "unfold" the com trajectory
        # we need an offset (multiples of the cell params) which needs to be incremented/decrementd
        # anytime the coordinate jumps by more then the cell param.
        offset = nm.zeros(com.shape, dtype="int32")
        current_offset = nm.zeros(com.shape[1:3], dtype="int32")
        self.pprint("Unfolding COM trajectory (non-parallel step!!!)")
        for i in xrange(1,self.nframes):
            d = com[i]-com[i-1]
            negative_jump = nm.less(d,-celldiag_half)
            positive_jump = nm.greater(d,celldiag_half)
            current_offset += nm.where(negative_jump, 1, 0)
            current_offset -= nm.where(positive_jump, 1, 0)
            offset[i,...] = current_offset
        self.pprint("done with unfolding")
        com_unfold = com + offset*celldiag
        # now lets do the msd
        msd = nm.zeros(com.shape, dtype="float64")
        entries = nm.zeros(com.shape[0], dtype="float64")
        entries[0] = 1.0
        if dt_max:
            # set the upper (unused) part of entries to one in order to avoid a division by zero
            entries[max_msd_frames+1:] = 1.0
        self.pprint("Now we compute the MSD up to a maximum deltat %10.5f ps" % (self.dt*max_msd_frames))
        for i in xrange(1+self.node_id,self.nframes,self.nodes):
            if (i%(self.nodes*100))==1: self.pprint("master node processing frame %d" % i)
            last_frame = i+max_msd_frames
            if last_frame > self.nframes: last_frame = self.nframes        
            r0 = com_unfold[i-1]
            rt = com_unfold[i:last_frame]
            dr = rt-r0
            msd[1:last_frame+1-i] += dr*dr
            entries[1:last_frame+1-i] += 1.0
        # now communicate and normalize
        if self.nodes>1:
            buf = nm.zeros(msd.shape, msd.dtype)
            self.comm.Allreduce(msd, buf, MPI.SUM)
            msd = buf
            buf = nm.zeros(entries.shape, entries.dtype)
            self.comm.Allreduce(entries, buf, MPI.SUM)
            entries = buf
        msd /= entries[:,nm.newaxis,nm.newaxis]
        if store:
            dataname = molname+"_MSD"
            self.pprint("storing MSD info in dataset %s" % dataname)
            if self.node_id == 0:
                self.msd_h5 = self.traj.require_dataset(dataname, shape=msd.shape, dtype=msd.dtype)
                self.msd_h5[...] = msd
            self.comm.Barrier()
        return msd

    def load_msd(self, msdname):
        self.msd_h5 = self.traj[msdname]
        return nm.array(self.msd_h5)
        
    def write_msd(self, msd, fname):
        self.pprint("\nWriting MSD data to %s" % fname)
        if self.node_id == 0:
            f = open(fname, "w")
            smsd = nm.sum(msd, axis=1)/float(msd.shape[1])
            tsmsd = nm.sum(smsd, axis=1)/3.0
            # convert from A^2 to nm^2
            smsd  *= 0.01
            tsmsd *= 0.01
            for i in xrange(self.nframes):
                f.write("%10.5f %12.6f %12.6f %12.6f %12.6f\n" % (i*self.dt*0.001, smsd[i,0], smsd[i,1], smsd[i,2], tsmsd[i]))
            f.close()
        self.comm.Barrier()
        return
        
    ####  oriantation stuff
    
    def compute_orient(self, molname, atompair, store=True):
        mollist = self.get_mollist(molname)
        nmols = len(mollist)
        self.pprint("\ncomputing the orientation vector of %d molecules with name %s" % (nmols,molname))
        self.pprint("   using atoms %d and %d to compute the vector" % (atompair[0], atompair[1]))
        self.pprint("Molecules are unwrapped .. this works only if the molecule is smaller then half the unit cell!!")
        atomlist = []
        for m in mollist:
            atomlist.append(self.mols[m])
        # allocate the array for the orientation
        # NOTE we could use the dtype of the xyz dataype here 
        orient = nm.zeros([self.nframes,nmols,3],dtype=self.xyz.dtype)
        if self.fixed_cell:
            celldiag = self.cell.diagonal()
            half_celldiag = celldiag/2.0
        for f in xrange(self.node_id,self.nframes, self.nodes):
            if (f%(self.nodes*100))==0: self.pprint("master node is processing frame %d" % f)
            if not self.fixed_cell:
                celldiag = self.cell[f].diagonal()
                half_celldiag = celldiag/2.0
            cxyz = self.xyz[f]
            for i,m in enumerate(atomlist):
                mxyz = cxyz[m]
                # unwrap the molecule, the reference image is the first atom:
                #         if any other atom is farther then the half cell we need to unwrap it
                #         NOTE this works only for molecules smaller then half the cell!!!
                dmxyz = mxyz[1:]-mxyz[0]
                mxyz[1:] += nm.where(nm.less   (dmxyz, -half_celldiag), 1.0, 0.0)*celldiag
                mxyz[1:] -= nm.where(nm.greater(dmxyz,  half_celldiag), 1.0, 0.0)*celldiag
                vect = mxyz[atompair[0]]-mxyz[atompair[1]]
                vect /= nm.sqrt(nm.sum(vect*vect))
                orient[f,i] = vect
        if self.nodes>1:
            buf = nm.zeros(orient.shape,dtype=orient.dtype)
            self.comm.Allreduce(orient, buf, MPI.SUM)
            orient = buf
        if store:
            dataname = molname+"_ORIENT"
            self.pprint("storing orientation info in dataset %s" % dataname)
            if self.node_id==0:
                self.orient_h5 = self.traj.require_dataset(dataname, shape=orient.shape, dtype=orient.dtype)
                self.orient_h5[...] = orient
            self.comm.Barrier()
        return orient
        
    def load_orient(self, molname):
        self.orient_h5 = self.traj[molname+"_ORIENT"]
        return nm.array(self.orient_h5)

    
        
    def orient_autocor(self, orient, molname, dt_max=None, store=True):
        """ computes the autocorrelation function of an orientation 
            assumes orient to be a (nframes, nmol, 3) array containing NORMALIZED orientation vectors"""
        self.pprint("\nComputing Orientational Autocorrelation Function (OACF) for molecule %s" % molname)
        # set up
        if dt_max:
            self.pprint("maximum deltat=%10.5f ps of a total of %10.5f ps" % (dt_max, self.nframes*self.dt))
            max_frames = int(dt_max/self.dt)
            if max_frames > self.nframes-1:
                raise ValueError, "requested delta t is more then the avialbale sampling time"
        else:
            max_frames = self.nframes
        nmols = orient.shape[1]
        oacf = nm.zeros([orient.shape[0],nmols], dtype="float64")
        entries = nm.zeros([orient.shape[0]], dtype="float64")
        entries[0] = 1.0
        if dt_max:
            # set the upper (unused) part of entries to one in order to avoid a division by zero
            entries[max_frames+1:] = 1.0
        for i in xrange(1+self.node_id,self.nframes, self.nodes):
            if (i%(self.nodes*100))==1: self.pprint("master node processing frame %d" % i)
            last_frame = i+max_frames
            if last_frame > self.nframes: last_frame = self.nframes        
            o0 = orient[i-1]
            ot = orient[i:last_frame]
            scalprod = nm.sum(ot*o0,axis=2)
            oacf[1:last_frame+1-i] += scalprod
            entries[1:last_frame+1-i] += 1.0
        # now normalize
        if self.nodes>1:
            buf = nm.zeros(oacf.shape, oacf.dtype)
            self.comm.Allreduce(oacf, buf, MPI.SUM)
            oacf = buf
            buf = nm.zeros(entries.shape, entries.dtype)
            self.comm.Allreduce(entries, buf, MPI.SUM)
            entries = buf
        oacf /= entries[:,nm.newaxis]
        oacf[0] = 1.0
        if store:
            dataname = molname+"_OACF"
            self.pprint("storing orientational autocorrelation function info in dataset %s" % dataname)
            if self.node_id==0:
                self.oacf_h5 = self.traj.require_dataset(dataname, shape=oacf.shape, dtype=oacf.dtype)
                self.oacf_h5[...] = oacf
            self.comm.Barrier()
        return oacf

    def load_oacf(self, molname):
        self.oacf_h5 = self.traj[molname+"_OACF"]
        return nm.array(self.oacf_h5)


    def write_oacf(self, oacf, fname, all=False, dt_max=None):
        self.pprint("\nWriting OACF to dat file")
        if self.node_id==0:
            nframes = oacf.shape[0]
            if dt_max:
                self.pprint("maximum deltat=%10.5f ps of a total of %10.5f ps" % (dt_max, nframes*self.dt))
                max_frames = int(dt_max/self.dt)
                if max_frames > nframes-1:
                    raise ValueError, "requested delta t is more then the avialbale sampling time"
            else:
                max_frames = nframes
            f = open(fname, "w")
            nmol = oacf.shape[1]
            soacf = nm.sum(oacf, axis=1)/float(nmol)
            for i in xrange(max_frames):
                line = "%10.5f " % (i*self.dt*0.001)
                if all:
                    line += (nmol*"%10.5f ") % tuple(oacf[i])
                line += "%12.6f \n" % soacf[i]
                f.write(line)
            f.close()
        return


    def map_com_2D(self, com, gridsize, name, store=True, orient=None):
        """ this is "cell based" in contrast to the "vertex based" code for the 3D map
            the reason is that VMD does it different for isosurface and volumeslice.
            please use this cube file ONLY for volume slices .. there is no z-info anyway"""
        self.pprint("\nComputing a 2D map of a com for molecule %s in the xy plane" % name)
        nframes = com.shape[0]
        grid = nm.zeros(gridsize, dtype="d")
        self.pprint("WARING: Map currently works only for constant volume ensembles!!")
        if not self.fixed_cell:
            raise ValueError, "map_com works currently only for fixed cell sizes"
        if not ((self.boundarycond==1)or(self.boundarycond==2)):
            raise ValueError, "Can't map for these boundaryconditions"
        # shift com by half a cell/2 up to get all coordinates above zeros
        celldiag = self.cell.diagonal()
        h = celldiag/(gridsize+[1])
        com += (celldiag)/2.0
        int_com = (com/h).astype("int32")
        entries = 0
        for f in xrange(self.node_id, nframes, self.nodes):
            for i, ind in enumerate(int_com[f].tolist()):
                if orient != None:
                    val = orient[f,i,2]
                    val *= val
                else:
                    val=1.0
                grid[tuple(ind[:2])]+=val
                entries += 1
        # comunicate 
        if self.nodes > 1:
            # first allreduce (also number of entries need to be allreduced)
            buf = nm.zeros(grid.shape, grid.dtype)
            self.comm.Allreduce(grid, buf, MPI.SUM)
            grid = buf
            entries = self.comm.allreduce(entries)
        # normalize and convert to a probability density
        grid /= (float(entries)*h.prod())
        if store:
            dataname = name+"_2DMAP"
            if self.nodes>1:
                self.pprint("storing 2Dmap in dataset %s" % dataname)
                self.map_h5 = self.traj.require_dataset(dataname, shape=grid.shape, dtype=grid.dtype)
                self.map_h5[...] = grid
        return grid

        
    def remove_from_pdlp(self, dataset):
        self.pprint("\nRemoving dataset %s from group /%s/traj" % (dataset, self.stagename))
        del(self.traj[dataset])
        return

