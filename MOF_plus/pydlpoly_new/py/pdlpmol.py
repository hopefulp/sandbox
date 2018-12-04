
"""
pdlpmol module implements the pdlpmol class for pydlpoly

part of the pydlpoly project
(C) 2014, CMC group, Rochus Schmid, Ruhr-Uni Bochum

This module and class works only in connection with pydlpoly and is meant to
keep the molecule methods seperate from the main code for better
maintainablility.

"""
import sys
import copy
import pydlpoly
import numpy
import numpy.random as nrand
import vectortools
from mpi4py import MPI

#dlp     = pydlpoly.dlp
#dlp_mol = pydlpoly.dlp_mol
#dlp_vdw = pydlpoly.dlp_vdw


# global flag for simplicitly .. this just mirrors the state of the lsplitvdw flag in the vdw module
split_vdw = False


class molecules:
    """
    This class holds the collection of all molecules in pydlpoly
    
    Just to keep methods applied to all molecules seperate from the rest of the code
    """
    def __init__(self, pd):
        """
        just get a reference to the pydlpoly class in order to get access tio everything else
        """
        self.pd = pd
        self.input_mol = self.pd.mol
        self.molnames = self.pd.mol.molnames
        #
        self.mlist = []
        for i in xrange(self.input_mol.nmols):
            self.mlist.append(mol(self.pd, i))
        # now generate the fortran datastructures
        self.nmolecules = self.input_mol.nmols  
        self.pd.dlp_mol.init_molecules(self.nmolecules)
        self.which = numpy.array(self.input_mol.whichmol,dtype="Int32")
        # add 1 to have fortran indexing!
        self.pd.dlp_mol.mol_which[:] = self.which+1
        mol_natoms = numpy.zeros([self.nmolecules], dtype="Int32")
        for i,m in enumerate(self.input_mol.mols):
            mol_natoms[i] = len(m)
        self.pd.dlp_mol.mol_natoms[:] = mol_natoms
        # use ekin per molecule by default
        self.pd.dlp_mol.lmol_ekin = True
        # if split vdw is requested in the pydlpoly class register here, too
        if self.pd.dlp_vdw.lspltvdw == True:
            split_vdw = True
            self.pd.pprint("Split lambda into dispersion and repulsion")
            self.pd.pprint("WARNING!!! This is an experimental feature and is not consistently implemented everywhere")
        return
    
    def get_list_from_name(self, molname):
        """ returns a list of molecule objects with the given name """
        mollist = []
        for m in self.mlist:
            if m.mname == molname: mollist.append(m) 
        return mollist
    
    def get_all_lambda_forces(self):
        """
        get a numpy array with all the forces
        """
        lamforce_vdw  = self.pd.dlp_mol.mol_dlam_vdw/self.pd.dlp.engunit
        lamforce_coul = self.pd.dlp_mol.mol_dlam_coul/self.pd.dlp.engunit
        return (lamforce_vdw, lamforce_coul)
        
    def set_all_lambda(self, mname, lvdw, lcoul):
        mollist = self.get_list_from_name(mname)
        self.pd.pprint("All molecules %s (%d) will be set to lambda %10.3f %10.3f" % (mname, len(mollist), lvdw, lcoul))
        for m in mollist: m.set_lambda(lvdw, lcoul)
        return
                
    def get_ekin(self):
        """gets the kinetic energy per molecule for all molecules
        
        .. todo::
        
            This should go into the molecuels class
        """
        ekin = self.pd.dlp_mol.mol_ekin/self.pd.dlp.engunit
        return ekin

    def use_lrcorr_lambda(self):
        """ switch on the use of lambda dependent long range vdw correction """
        self.pd.dlp_mol.lmol_lrcorr = True
        return
        
    def get_lrcorr_lambda_force(self):
        """ gets the force on lambda from the long range correction
        
        WARNING: this works ONLY if a sinlge molecule has a lambda between 0 and 1
        (balck sheep) the force is computed under this assumption """
        lamforce_vdw_lrcorr = self.pd.dlp_mol.mol_dlam_lrcorr/self.pd.dlp.engunit
        return lamforce_vdw_lrcorr
        
        

class mol:
    """
    This class provides molecule level methods for pydlpoly
    
    Each pdlpmol object is a molecule within pydlpoly. They are automatically generated upon
    setup in pydlpoly and a number of things have to stay in pydlpoly.
    This class is really not meant to be used by itself but always in connection with pydlpoly.
    The main purpose is to keep these methods seperate from pydlpoly.py
    
    .. note::
    
        For Developers:
        The information about molecules (names/atoms belongign to a molecule etc)
        is generated and maintained in tha :class:assign_FF.mol in the `assing_FF.py` module.
        For historic reasons it will stay there. Also IO to/from pdlpio (for example in a restart)
        will put the data there. The objects here also contain the same info and are just to have a simple
        handle to do some operations on these molecules.  
    
    """
    
    def __init__(self, pd, molid):
        """
        Generates a pydlpoly molecule object
        
        A molecule is a system wich is not bonded to any other atom. This means a molecule interacts with all other
        atoms *only* by non-bonded interactions (Coulombic and vdW)
        
        In contrast to the previous implementation, where the molecule info was kept in the
        :class:assign_FF.mol class, all info is pulled over into this class, but the objects
        are still generated py pydlpoly to set up the internal F90 data structures.
        
        There are two main reasons to use molecule objects:
            - Move a complete molecule (translate, rotate)
            - change its lambda values (vdw and Coulomb) to switch its interaction on and off 
              (also the kinetic energy can be requested or velocities can be initialized or quenched)
        
        :Parameters:
            - pd :   pydlpoly class to access all necesary information
            - molid: the molid (integer)
        
        .. note::
        
            For operations like translate_com and rotate it is necessary (and it is simply assumed)
            that the atoms of the molecule are packed (all in same periodic image).
            We need a method to enforce this when using these methods after an MD or a minimization where 
            the atoms can potentially be wrapped (and unpacked) by the image subroutines in pydlpoly
        
        """
        
        self.pd      = pd
        self.molid   = molid
        self.matoms  = self.pd.mol.mols[self.molid]
        self.natoms  = len(self.matoms)
        self.mtype   = self.pd.mol.moltypes[self.molid]
        self.mname   = self.pd.mol.molnames[self.mtype]
        self.mdegfre = len(self.matoms)*3
        self.masses  = self.pd.get_subset_masses(self.matoms)
        self.tot_mass= self.masses.sum()
        # charge for use in QEq
        self.Q       = 0.0
        # switch to splitvdw versions
        if self.pd.dlp_vdw.lspltvdw == True:
            self.get_lambda_force = self.get_lambda_force_splitvdw
            self.set_lambda       = self.set_lambda_splitvdw
        return
    

    def __str__(self):
        return self.mname
    
    def get_lambda_force(self):
        """
        Return the lambda force for this molecule
        
        :Returns:
            dlambda vdW, dlambda_Coulomb
            
        """
        dlmb_vdw = self.pd.dlp_mol.mol_dlam_vdw[self.molid]/self.pd.dlp.engunit
        dlmb_coul= self.pd.dlp_mol.mol_dlam_coul[self.molid]/self.pd.dlp.engunit
        return (dlmb_vdw, dlmb_coul)
        
    def get_lambda_force_splitvdw(self):
        """
        Return the lambda force for this molecule
        
        :Returns:
            dlambda vdW, dlambda_Coulomb, dlambda vdW repulsion
            
        """
        dlmb_vdw = self.pd.dlp_mol.mol_dlam_vdw[self.molid]/self.pd.dlp.engunit
        dlmb_coul= self.pd.dlp_mol.mol_dlam_coul[self.molid]/self.pd.dlp.engunit
        dlmb_vdwr = self.pd.dlp_mol.mol_dlam_vdwr[self.molid]/self.pd.dlp.engunit
        return (dlmb_vdw, dlmb_coul, dlmb_vdwr)


    def set_lambda(self, vdw, coul):
        """
        Set the lambda of the molecule
        
        :Parameters:
            - vdw  [float 0.0-1.0]
            - coul [float 0.0-1.0]
            
        """
        self.pd.dlp_mol.mol_lamb_vdw[self.molid] =vdw
        self.pd.dlp_mol.mol_lamb_coul[self.molid]=coul
        # set flag to enforce recalc of energy
        self.pd.set_atoms_moved()
        return
        
    def set_lambda_splitvdw(self, vdw, coul, vdwr=0.0):
        """
        Set the lambda of the molecule
        
        :Parameters:
            - vdw  [float 0.0-1.0]
            - coul [float 0.0-1.0]
            - vdwr  [float 0.0-1.0]
            
        """
        self.pd.dlp_mol.mol_lamb_vdw[self.molid] =vdw
        self.pd.dlp_mol.mol_lamb_coul[self.molid]=coul
        self.pd.dlp_mol.mol_lamb_vdwr[self.molid] = vdwr
        # set flag to enforce recalc of energy
        self.pd.set_atoms_moved()
        return


    def set_lambda_ekin(self, lambekin):
        """
        Set molecules lambda for ekin
        """
        self.pd.dlp_mol.mol_lamb_ekin[self.molid] = lambekin
        return

    def get_xyz(self):
        return self.pd.get_subset_xyz(self.matoms)
    
    def get_vel(self):
        return self.pd.get_subset_vel(self.matoms)
    
    def set_xyz(self, xyz):
        self.pd.set_subset_xyz(xyz, self.matoms)
        self.pd.rebuild_nlist()
        return
        
    def get_force(self):
        return self.pd.get_subset_force(self.matoms)
        
    def add_force(self, force):
        self.pd.add_subset_force(force, self.matoms)
        return
                                
    def get_com(self):
        """ compute the center of mass for a molecule
        NOTE: the molecule is wrapped into the box first using the first atom as the reference image """
        bcond = self.pd.get_bcond()
        if (bcond > 0) and (bcond < 3):
            # the system is orthormobic at most -> we can live with the celldiagonal
            celldiag = self.pd.get_cell().diagonal()
            xyz = self.get_xyz()
            dxyz = xyz[1:]-xyz[0]
            xyz[1:] += numpy.where(numpy.less_equal(dxyz, -celldiag*0.5), 1.0, 0.0)*celldiag
            xyz[1:] -= numpy.where(numpy.greater   (dxyz,  celldiag*0.5), 1.0, 0.0)*celldiag
            com = numpy.sum(xyz*self.masses[:,numpy.newaxis], axis=0)/self.tot_mass
            return com
        elif bcond == 3:
            # triclinic
            cell = self.pd.get_cell()
            icell = numpy.linalg.inv(cell)
            pos = self.get_xyz()
            fpos = numpy.dot(pos, icell)
            dfpos = fpos[1:]-fpos[0]
            fpos[1:] += numpy.where(numpy.less_equal(dfpos, -0.5), 1.0, 0.0)
            fpos[1:] -= numpy.where(numpy.greater   (dfpos,  0.5), 1.0, 0.0)
            pos = numpy.dot(fpos, cell)
            com = numpy.sum(pos*self.masses[:,numpy.newaxis], axis=0)/self.tot_mass
        else:
            raise ValueError, "get_com not implemented for this bcond"
        return
        
    def translate_mol(self, shift):
        """shifts a molecule by vector
        check if wrapping works ok ... also we might have to enforce neigbor list update"""
        xyz = self.get_xyz()
        xyz += shift
        self.set_xyz(xyz)
        return

    #def translate_PBC(self, m, shift):
        ## shifts a molecule by vector taking to account periodic boundary conditions
        #celldiag = self.get_cell().diagonal()
        #matoms = self.mol.mols[m]
        #xyz = self.get_subset_xyz(matoms)
        #xyz += shift
        #center = numpy.array([0.0,0.0,0.0])
        #dxyz = xyz[0:]-center
        #xyz[0:] += numpy.where(numpy.less(dxyz, -celldiag*0.5), 1.0, 0.0)*celldiag
        #xyz[0:] -= numpy.where(numpy.greater   (dxyz,  celldiag*0.5), 1.0, 0.0)*celldiag
        #self.set_subset_xyz(xyz, matoms)
        #self.rebuild_nlist()
        #return
    
    def move_com_to(self, pos):
        """ moves molecule COM to a given position  """
        self.translate_mol(pos-self.get_com())
        return


    def rotate_random(self):
        xyz = self.get_xyz()
        # rotate molecule by random quaternion
        xyz = vectortools.rotate_random(xyz)
        #celldiag = self.get_cell().diagonal()
        #center = numpy.array([0.0,0.0,0.0])
        #dxyz = xyz[0:]-center
        #xyz[0:] += numpy.where(numpy.less(dxyz, -celldiag*0.5), 1.0, 0.0)*celldiag
        #xyz[0:] -= numpy.where(numpy.greater   (dxyz,  celldiag*0.5), 1.0, 0.0)*celldiag
        self.set_xyz(xyz)
        return

    
    def quench(self, outf=None):
        """ quenches the velocities of this molecule (and removes it from the degrees of freedom for the thermostat)
        NOTE: currently this does not touch the COM velocities of the total system
              in other words doing this will lead to a COM (rot)momentum of the total system"""
        vel = self.pd.get_subset_vel(self.matoms)
        ekin = 0.5*numpy.sum(self.masses[:,numpy.newaxis]*vel*vel)/self.pd.dlp.engunit
        vel[::] = 0.0
        self.pd.set_subset_vel(vel, self.matoms)
        #self.pd.set_thermostat_sigma(self.pd.get_degfree()-float(self.mdegfre))
        self.pd.set_degfree(self.pd.get_degfree()-float(self.mdegfre))
        self.set_lambda_ekin(0.0)
        self.pprint("  QUENCHed molecule %3d, %3d degrees of freedom with %10.5f kcal/mol kinetic energy removed" % \
             (self.molid, self.mdegfre, ekin), outf=outf)
        self.pprint("     remaining degrees of freedom: %d " % int(self.pd.get_degfree()), outf=outf)
        return ekin
        
    def init_vel(self, exact=True, outf=None):
        """ this method reinitalizes the velocities of the given molecule to (and adds it back to the degrees of freedom for the thermostat)
        NOTE: the COM velocity is not removed ... nonsense for an atom anyway 
            this needs to be tested """
        T      = self.pd.get_temperature()
        new_vel = nrand.normal(0.0, numpy.sqrt(pydlpoly.boltz*T/self.masses[:,numpy.newaxis]), size=[self.natoms,3])
        if self.pd.nodes > 1:
            # we are parallel and communicate the random from node 0
            self.pd.local_comm.Bcast([new_vel, MPI.DOUBLE], root=0)
        if exact:
            # rescale to the exact temp
            init_Ekin = 0.5*numpy.sum(self.masses[:,numpy.newaxis]*new_vel*new_vel)
            new_vel *= numpy.sqrt(0.5*3*self.natoms*pydlpoly.boltz*T/init_Ekin)
        self.pd.set_subset_vel(new_vel, self.matoms)
        ekin = 0.5*numpy.sum(self.masses[:,numpy.newaxis]*new_vel*new_vel)/self.pd.dlp.engunit
        #self.pd.set_thermostat_sigma(self.pd.get_degfree()+float(self.mdegfre))
        self.pd.set_degfree(self.pd.get_degfree()+float(self.mdegfre))
        self.set_lambda_ekin(1.0)
        self.pprint("  INIT velocities of molecule %3d with a kinetic energy of %10.5f kcal/mol" % \
            (self.molid, ekin), outf=outf)
        return ekin

    def set_charge(self, Q):
        self.Q = Q
        return
                           
    def pprint(self, out, outf=None):
        if outf:
            outf.write(out+"\n")
        else:
            self.pd.pprint(out)  
        return
        
    