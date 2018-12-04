# -*- coding: utf-8 -*-

import numpy
from molsys.util import elems as elements
import string
import copy

class hessian(object):

    def __init__(self, pd, lmp = False):
        self.pd = pd
        self.lmp = lmp
        if lmp:
            self.get_force = self.calc_force
            self.pd.lmps.command('run 0 pre yes post no')
        else:
            self.get_force = self.pd.calc_energy_force
        return

    def calc_force(self):
        self.pd.lmps.command('run 1 pre no post no')
        return None, self.pd.get_force()

class doublehessian(hessian):

    def __call__(self, disp=0.001):
#        self.pd.pprint("calculation of MM hessian started")
        if self.lmp: self.pd.lmps.command('run 0 pre yes post no')
        delta = disp
        n = self.pd.get_natoms()
        hessian=numpy.zeros((3*n,3*n), dtype="float64") 
        xyz_orig = self.pd.get_xyz()
        xyz = xyz_orig.copy()
        if self.pd.QMMM:
            # if this is a QMMM run then we skip the QM calc in calc_energy()
            self.pd.QMMM_interface.skip_QM_energy = True
        for i in range(n): 
            for j in range(3):
                ### positiver gradient ###
                xyz[i,j] = xyz_orig[i,j] + delta      
                self.pd.set_xyz(xyz)
                #self.set_atoms_moved()
                e,fp = self.get_force()
                ### negativer gradient ###
                xyz[i,j] = xyz_orig[i,j] - delta
                self.pd.set_xyz(xyz)
                #self.set_atoms_moved()
                e,fm = self.get_force()
                xyz[i,j] = xyz_orig[i,j]
                ### insgesamter gradient ###
                fi = (fm - fp)/(2*delta)
                hessian[:,3*i+j]=fi.ravel()
#        self.pd.pprint("finite difference hessian done")
        if self.pd.QMMM:
            # now switch QM energy back on
            self.pd.pprint("calculation of QM hessian started")
            self.pd.QMMM_interface.skip_QM_energy = False
            # now do analytic QM-Hessian
            hessian += self.pd.QMMM_interface.build_transhessian()
            self.pd.pprint("calculation of QM hessian done")
        self.s_hessian = (hessian + hessian.transpose())/2 ### kcal/(A**2 * mol) ###
        self.pd.set_xyz(xyz_orig)
        if self.lmp: self.pd.lmps.command('run 0 pre yes post no')
        return self.s_hessian

class singlehessian(hessian):

    def __call__(self, disp=0.001):
#        self.pd.pprint("calculation of MM hessian started")
        if self.lmp: self.pd.lmps.command('run 0 pre yes post no')
        delta = disp
        n = self.pd.get_natoms()
        hessian=numpy.zeros((3*n,3*n), dtype="float64") 
        xyz_orig = self.pd.get_xyz()
        e, f = self.get_force()
        forig = f.ravel()
        xyz = xyz_orig.copy()
        for i in range(n): 
            for j in range(3):
                ### positiver gradient ###
                xyz[i,j] = xyz_orig[i,j] + delta      
                self.pd.set_xyz(xyz)
                e,fp = self.get_force()
                xyz[i,j] = xyz_orig[i,j]
                ### insgesamter gradient ###
#                fi = - fp/delta
                hessian[:,3*i+j]=(forig-fp.ravel())/delta
#        self.pd.pprint("finite difference hessian done")
        self.s_hessian = (hessian + hessian.transpose())/2 ### kcal/(A**2 * mol) ###
        self.pd.set_xyz(xyz_orig)
        if self.lmp: self.pd.lmps.command('run 0 pre yes post no')
        return self.s_hessian

class hessutils(object):

    def __init__(self, xyz, hessian, elems, masses = None):
        self.s_hessian = hessian
        self.xyz = xyz
        if len(elems) != numpy.shape(self.s_hessian)[0]/3:
            print 'Shapes does not match!'
            raise IOError
        self.natoms = len(elems)
        self.elems = []
        if masses is not None: 
            self.masses = masses
        else:
            self.masses = []
        for i in elems:
            e = i.lower().split()[0]
            self.elems.append(e)
            if masses is None:
                self.masses.append(string.atof(elements.mass[e]))

    @staticmethod
    def get_massmatrix(natoms, masses): 
        m = numpy.zeros((natoms,3), dtype="float64")
        m[:,0] = m[:,1] = m[:,2] = 1/(numpy.sqrt(masses))
        mm = numpy.outer(m.ravel(),m.ravel())
        return mm*4.184e26


    def calc_mass_weighted_hessian(self):
        m = numpy.zeros((self.natoms,3), dtype="float64")
        m[:,0] = m[:,1] = m[:,2] = 1/(numpy.sqrt(self.masses))
        mm = numpy.outer(m.ravel(),m.ravel())
        mw_hessian = mm * self.s_hessian * 4.184e26 ###scale factor from kcal/(A**2 * g) to SI (1/s**2) ###
        return mw_hessian
    
    def calc_eigenvalues(self):
        c=3.0*10.0**8.0
        freq = numpy.zeros((3*self.natoms), dtype="float64")
        evalue, evector = numpy.linalg.eigh(self.calc_mass_weighted_hessian()) 
        for i in xrange(3*self.natoms):
            if evalue[i] >= 0:
                freq[i] = numpy.sqrt(evalue[i])/(2*numpy.pi)
            else:
                freq[i] = (-1) * numpy.sqrt(abs(evalue[i]))/(2*numpy.pi)
        l = c/freq
        wn = ((1/l)/100)
        self.wn = wn
        self.evector = evector
        return wn, evector

    def screw(self, vector, factor):
        xyz = copy.deepcopy(self.xyz)
        mass = 1.0/(numpy.sqrt(self.masses))
        return xyz
    
    def write_molden_input(self, fname):
        xyz = self.xyz * (1.0/0.52917721092)
        freq, evector = self.calc_eigenvalues()
        evector = evector * (1.0/0.52917721092) 
        mass = 1/(numpy.sqrt(self.masses)) 
        minput = open(fname, "w")
        ### Title section ###
        minput.write("[Molden Format]\n[Title]\n\n")
        ### coord section ###
        minput.write("[Atoms] Au\n")
        a = 1
        for i in xrange(self.natoms):
            minput.write("%s %s %s %12.6f %12.6f %12.6f\n" % (self.elems[i], a, elements.number[self.elems[i]], xyz[i,0], xyz[i,1], xyz[i,2]))
            a += 1
        ### freq section ###
        minput.write("[FREQ]\n")
        for j in range(numpy.shape(freq)[0]):
            minput.write("%12.6f\n" % (freq[j]))
        ### fr-coord section 
        minput.write("[FR-COORD]\n")
        for k in xrange(self.natoms):
            minput.write("%s %12.6f %12.6f %12.6f\n" % (self.elems[k], xyz[k,0], xyz[k,1], xyz[k,2])) 
        ### fr-norm-coord section ###
        minput.write("[FR-NORM-COORD]\n")
        for l in range(numpy.shape(freq)[0]):
            minput.write("vibration %s\n" % (l+1)) # Normal modes
            for m in range(len(mass)):
                minput.write("%12.6f %12.6f %12.6f\n" % (evector[3*m+0,l]*mass[m], evector[3*m+1,l]*mass[m], evector[3*m+2,l]*mass[m]))
        minput.close()
        return




