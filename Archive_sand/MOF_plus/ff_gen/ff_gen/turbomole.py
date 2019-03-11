#! /usr/bin/env python
########################################
#
#  pysthon script to control turbomole
#
########################################

import numpy
import os
import subprocess
import shutil
import string
import time
import elements as elems

BOHR2A = 0.5291772
HA2KCAL = 627.5095

class turbomole:

#    def __init__(self, xyz, elements, sed = None, protfile = None):
    def __init__(self, natoms, elements=None, sed = None, protfile = None):
        self.sed = sed
        self.protfile = protfile
        self.elements = elements
#        self.xyz = xyz
        self.natoms = natoms
#        self.masses = []
#        for i in range(self.natoms):
#            self.masses.append(string.atof(elems.mass[string.split(self.elements[i].lower())[0]]))
        return

    def write_coord(self, xyz):
        BOHR2A  = 0.5291772
        f = open("coord", "w")
        f.write("$coord\n")
        for i in range(self.natoms):
            c = xyz[i]/BOHR2A
            e = string.lower(self.elements[i])
            f.write("%15.10f %15.10f %10.10f %s\n" % (c[0], c[1], c[2], e))
        f.write("$end\n")
        f.close()
        return

    def setup(self):
        print "define is running"
        os.system ("define < %s > define.out" % self.protfile)
        time.sleep(2.0)
        if self.sed:
            os.system("sed '%s' control > control.sed" % self.sed)
            shutil.copy("control.sed", "control")
        return

    def optimize(self):
        ferr = open("jobex_blurp.dump","w")
        f = open("jobex.out","w")
        print "jobex is running"
        retcode = subprocess.call(["jobex", "-ri"], stdout=f, stderr=ferr)
        if retcode != 0: raise IOError, "Error: jobex failed!!!"
        f.close()
        ferr.close()
        print "jobex has finished"
        os.sytem('t2xyz|tail -%s > final.xyz' % self.natoms +2)
        return

    def calc_energy(self, force):
        self.cycle += 1
        if os.path.exists("control"):
            #self.cycle += 1
            ferr = open("turbo_blurb.dump","w")
            f = open("turbo_%d.out" % self.cycle, "w")
            retcode = subprocess.call(["ridft"], stdout=f, stderr=ferr)
            if retcode != 0: raise IOError, "QMMM Error: RIDFT failed!!!"
            if force:
                retcode = subprocess.call(["rdgrad"], stdout=f, stderr=ferr)
                if retcode != 0: raise IOError, "QMMM Error: RDGRAD failed!!!"
            f.close()
            ferr.close()
            if force:
                f = self.read_gradient()
            else:
                e = self.read_energy()
        return e, f

    def calc_hessian(self):
        print("preparing the control file")
        os.system("sed -e '/^$ricore/d' -e '/^$rij/d' -e '/^$jbas/d' -e 's/$marij/$noproj/g'  control > control.aoforce")
        shutil.copy("control.aoforce", "control")
        ferr = open("aoforce_blurp.dump","w")
        f = open("aoforce.out","w")
        print("aoforce is running")
        retcode = subprocess.call(["aoforce"], stdout=f, stderr=ferr)
        if retcode != 0: raise IOError, "Error: aoforce failed!!!"
        f.close()
        ferr.close()
        print("aoforce done")
        os.system('tm2molden < input')
        return

    def read_optimized_xyz(self):    
        os.system('t2xyz|tail -%s > final.xyz' % (self.natoms +2))
        f = open("final.xyz", 'r')
        line = string.split(f.readline())
        line = string.split(f.readline())
#        if int(line[0]) != self.natoms:
#            raise IOError, "Problem with the final geometry, number of atoms does not match!"
        xyz = numpy.zeros((self.natoms, 3))
        for i in range(self.natoms):
            line = string.split(f.readline())
            for j in range(3):
                xyz[i,j] = line[j+1]
        self.xyz = xyz
        print 'Optimized geometry was read'
        return

    def read_gradient(self, tcycle):
        f = open("gradient", "r")
        line = string.split(f.readline())
        if line[0] != "$grad":
            raise IOError, "Not a TURBOMOLE gradient file!"
        line = string.split(f.readline())
        cycle = string.atoi(line[2])
        while cycle < tcycle:
            for i in xrange(2*(self.natoms)):
                f.readline()
            line = string.split(f.readline())
            cycle = string.atoi(line[2])
        # now we are at the right spot .. get the energy first
        energy = string.atof(line[6])*HA2KCAL
        # skip the coords
        for i in xrange(self.natoms):
            f.readline()
        gradient = []
        for i in xrange(self.natoms):
            rawline = f.readline().replace("D", "E")
            line = string.split(rawline)
            gradient.append(map(string.atof, line[:3]))
        gradient = numpy.array(gradient)
        gradient *= HA2KCAL/BOHR2A
        return energy, gradient
            
    def read_energy(self):
        f = open("energy", "r")
        line = string.split(f.readline())
        if line[0] != "$energy":
            raise IOError, "Not a TURBOMOLE energy file!"
        line = string.split(f.readline())
        cycle = string.atoi(line[0])
        while cycle < self.cycle:
            line = string.split(f.readline())
            cycle = string.atoi(line[0])
        # now we are at the right spot .. get the energy
        energy = string.atof(line[1])*HA2KCAL
        return energy

    def read_hessian(self, fname = "control", keyword ="$nprhessian"):
        n = self.natoms
        f = open(fname, "r")
        found = False
        while not found:
            line = string.split(f.readline())
            try:
                if line[0] == keyword: found = True
            except IndexError:
                print 'No hessian found in the control file using the keyword: ' ,keyword
                raise IndexError
        hessian = numpy.zeros((3*n,3*n), dtype = "float64")
        line = string.split(f.readline())
        r = int(round(((3*n/5.0 - 3*n/5)*5)))
        if r == 0:
            a = 0
        else:
            a = 1
        for x in xrange(3*n):
            for y in xrange((3*n/5) + a): ### a reinbasteln
                l = len(line)
                if y+1 <= 3*n/5:
                    values = numpy.array(map(string.atof, line[l-5:]))
                    hessian[x,(y*5):(y*5)+5] = values
                else:
                    values = numpy.array(map(string.atof, line[l-r:]))
                    hessian[x,(y*5):(y*5)+r] = values ###????
                line = string.split(f.readline())
        f.close()
        self.hessian = (hessian* 2240.877106) ### scale factor from a.u. to kcal/(A**2 * mol)  
        return

    def read_kollman(self,fname = "charge.out"):
        f = open(fname, "r")
        found = False
        while not found:
            line = string.split(f.readline())
            try:
                if line == ["charges","resulting","from","fit:"]: found = True
            except IndexError:
                print 'No charges found in the file: ' , fname
                exit()
        for i in range(3): line = f.readline()
        if line.split()[0] != 'atom': f.readline()
        charges = []
        for i in range(self.natoms):
            line = string.split(f.readline())
            charges.append(float(line[3]))
        self.charges = numpy.array(charges)
        f.close()
        return

    def read_control(self, fname, key = "$last step"):
        control = ""
        keylen = len(key)
        with open(fname, "r") as f:
            for line in f:
                if line[:keylen] != key:
                    control += line
                else:
                    break
        return control





