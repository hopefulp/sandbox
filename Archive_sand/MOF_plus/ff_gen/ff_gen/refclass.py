#!/usr/bin/env python
# -*- coding: utf-8 -*-

###################################################
#
#  Reference class for storing the reference data
#
###################################################

import numpy
import os
import copy
import turbomole
import h5py
import string
import elements as elems
import pdlpio
import shutil

class refclass2(object):

    def __init__(self, fname):
        self.f = h5py.File(fname, 'r')
        self.name = string.join(fname.split('/')[-1].split('.')[:-1],'.')

    def __call__(self, info, branch = 'hessians', tag = 'primary'):
        if branch == 'system':
            return self.f[branch][info].value
        else:
            return self.f[branch][tag][info].value



class refclass:

    def __init__(self, refname, refdata, refprog = 'turbomole', protfile = None, sed = None, create_ref=False):
        self.refname = refname
        self.xyz = refdata[0]
        self.elements = refdata[1]
        self.natoms = numpy.shape(self.xyz)[0]
        self.protfile = protfile
        self.sed = sed
        if os.path.isfile('%s.hdf5' % self.refname):
            self.f = h5py.File('%s.hdf5' % self.refname, 'r+')
            self.natoms = self.__call__('natoms', 'system')
            self.elements = self.__call__('elements', 'system')
            self.exist = True
        else:
            self.f = h5py.File('%s.hdf5' % self.refname, 'a')
            self.exist = False 
            self.create_reference()
            if create_ref == True:
                self.read_hessian('primary')
        return


    def create_reference(self):
#        retval = os.getcwd()
#        os.chdir('ref')
        self.f.create_group('system')
        self.f['system'].create_dataset('natoms', data = self.natoms)
        self.f['system'].create_dataset('elements', data = self.elements)
#        os.chdir(retval)
        return


    def read_structures(self, nstruc, tag = 'potential'):
        retval = os.getcwd()
        os.chdir('ref/strucs')
        self.nstruc = nstruc
        self.f.require_group('forcematch')
        self.f['forcematch'].require_group(tag)
        self.f['forcematch'][tag].require_dataset('structures',
                shape = (self.nstruc, self.natoms, 3), dtype = 'float64')
        self.f['forcematch'][tag].require_dataset('forces',
                shape = (self.nstruc, self.natoms, 3), dtype = 'float64')
        self.f['forcematch'][tag].require_dataset('energies',
                shape = (self.nstruc,1), dtype = 'float64')
        for i in range(nstruc):
            self.f['forcematch'][tag]['structures'][i,:,:]= self.read_xyz('%s.xyz' % i)
        os.chdir(retval)            
        return

    def allocate_structures(self, nstruc, tag = 'md', vectors = False):
        self.f.require_group('forcematch')
        self.f['forcematch'].require_group(tag)
        self.f['forcematch'][tag].require_dataset('structures',
                shape = (nstruc, self.natoms, 3), dtype = 'float64')
        self.f['forcematch'][tag].require_dataset('forces',
                shape = (nstruc, self.natoms, 3), dtype = 'float64')
        self.f['forcematch'][tag].require_dataset('energies',
                shape = (nstruc,1), dtype = 'float64')
        if vectors:
            self.f['forcematch'][tag].require_dataset('vectors',
                    shape = (nstruc, self.natoms, 3), dtype = 'float64')
        return

    def read_structures_from_pdlp(self, fpdlp, tag = 'md', dir='ref/strucs'):
        retval = os.getcwd()
        os.chdir(dir)
        pdio = pdlpio.pdlpio(fpdlp)
        if pdio.has_stage('ref') == False:
            print 'No reference stage in %s' % fpdlp
            raise IOError
        if pdio.h5file['ref'].keys().count('traj') == 0:
            print 'No trajectory in reference stage in %s' % fpdlp
            raise IOError
        self.nstruc = numpy.shape(pdio.h5file['ref']['traj']['xyz'].value)[0]
        self.f.require_group('forcematch')
        self.f['forcematch'].require_group(tag)
        self.f['forcematch'][tag].require_dataset('structures',
                shape = (self.nstruc, self.natoms, 3), dtype = 'float64')
        self.f['forcematch'][tag].require_dataset('forces',
                shape = (self.nstruc, self.natoms, 3), dtype = 'float64')
        self.f['forcematch'][tag].require_dataset('energies',
                shape = (self.nstruc,1), dtype = 'float64')
        self.f['forcematch'][tag]['structures'][:,:,:] = copy.deepcopy(
                pdio.h5file['ref']['traj']['xyz'].value)
        os.chdir(retval)
        return

    def read_hessian(self, tag = 'secondary', fname = 'control', keyword = '$nprhessian',
            fcharges = None, path = None):
        retval = os.getcwd()
        if path == None:
            os.chdir('ref/%s' % tag)
        else:
            os.chdir(path)
        tb = turbomole.turbomole(self.natoms, self.elements)
        tb.read_optimized_xyz()
        tb.read_hessian(fname, keyword)
        if fcharges is not None:
            tb.read_kollman(fcharges)
            charges = tb.charges
            print 'Charges have been read'
        else:
            charges = None
        print 'Hessian has been read'
        os.chdir(retval)
        self.write_hessian(tag, tb.xyz, tb.hessian, charges)
        return


    def write_hessian(self, tag, coord, hessian, charges = None):
        assert numpy.shape(coord)[0] == self.natoms
        assert numpy.shape(coord)[1] == 3
        assert numpy.shape(hessian)[0] == 3*self.natoms
        assert numpy.shape(hessian)[1] == 3*self.natoms
        self.f.require_group('hessians')
        self.f['hessians'].require_group(tag)
        self.f['hessians'][tag].require_dataset('energy', shape = (1,1), dtype = 'float64')
        self.f['hessians'][tag].require_dataset('coord', shape = (self.natoms, 3), dtype = 'float64')
        self.f['hessians'][tag].require_dataset('hessian', shape = (3*self.natoms, 3*self.natoms), dtype = 'float64')
        self.f['hessians'][tag]['coord'][:] = coord
        self.f['hessians'][tag]['hessian'][:] = hessian
        if charges is not None:
            self.add_charges(tag,charges)
        return

    def write_grid(self, tag, ngrid, grid, energies):
        self.f.require_group('grids')
        self.f['grids'].require_group(tag)
        self.f['grids'][tag].require_dataset('grid', shape = (ngrid, 3), dtype = 'float64')
        self.f['grids'][tag].require_dataset('energy', shape = (1,ngrid), dtype = 'float64')
        self.f['grids'][tag].require_dataset('ngrid', shape = (1,1), dtype = 'float64')
        self.f['grids'][tag]['grid'][:] = grid
        self.f['grids'][tag]['ngrid'][:] = ngrid
        self.f['grids'][tag]['energy'][:] = energies
        return

    def add_charges(self,tag,charges):
        assert numpy.shape(charges)[0] == self.natoms
        self.f['hessians'][tag].require_dataset('charges', shape = (1, self.natoms), dtype = 'float64')
        self.f['hessians'][tag]['charges'][:] = charges
        return

    def calculate_structures(self, tag, startstruc = None, endstruc = None, restart = False):
        tb = turbomole.turbomole(self.natoms, self.elements, protfile='../../../'+self.protfile, sed=self.sed)
        retval = os.getcwd()
        if startstruc == None:
            startstruc = 0
        if endstruc == None:
            endstruc = numpy.shape(self.f['forcematch'][tag]['structures'].value)[0]
        if restart == True:
            if os.path.isdir('ref/strucs/restart') != True:
                print 'Restart directory not found'
                raise IOError
        for i in range(startstruc, endstruc):
            if restart == True:
                shutil.copytree('ref/strucs/restart', 'ref/strucs/%s' %i)
                os.chdir('ref/strucs/%s' %i)
                tb.write_coord(self.f['forcematch'][tag]['structures'][i,:,:])
            else:
                os.mkdir('ref/strucs/%s' % i)
                os.chdir('ref/strucs/%s' % i)
                tb.write_coord(self.f['forcematch'][tag]['structures'][i,:,:])
                tb.setup()
            self.write_jobscript(i, tag, retval)
            os.system('qsub q_%s' % i)
            os.chdir(retval)
        return

    def write_jobscript(self, i, tag,  path):
        f = open("q_%s" % i, "w")
        f.write("#!/bin/bash\n\n")
        f.write("export PARNODES=2\n")
        f.write("module load ffgen_intel")
        f.write("#$ -S /bin/bash\n")
        f.write("#$ -cwd\n")
        f.write("#$ -q par.q\n")
        f.write("#$ -pe shmem 2\n")
        f.write("#$ -j y\n")
        f.write("#$ -V\n")
        f.write("#$ -N %s_%s\n" % (tag, i))
        f.write("export struc=%s\n" % i)
        f.write("export tag=%s\n" % tag)
        f.write("ridft > ridft.out\n")
        f.write("rdgrad > grad.out\n")
        hpath = path + '/%s.hdf5' % self.refname
        f.write("collect_data.py -i %s > collect.out" % hpath) 
        return

    def read_xyz(self,fname):
        f = open(fname, 'r')
        line = string.split(f.readline())
        line = string.split(f.readline())
        xyz = numpy.zeros((self.natoms, 3))
        for i in range(self.natoms):
            line = string.split(f.readline())
            for j in range(3):
                xyz[i,j] = line[j+1]
        return xyz


    def __call__(self, info, branch = 'hessians', tag = 'primary'):
        if branch == 'system':
            return self.f[branch][info].value
        else:
            return self.f[branch][tag][info].value
