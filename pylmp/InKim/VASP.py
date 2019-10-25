#!/home/noische/python

import sys
import os
import getopt
import math
import time
import re

import numpy as np
from numpy import sin, cos, sqrt, arccos, fabs, pi, rad2deg, deg2rad
import bgf
import bgftools
import nutils as nu

usage = """This is a VASP class in python.
You can use this class as follows:
>> from VASP import *

- POSCAR file import
>> myfile = VASP("POSCAR")

- Save POSCAR as BGF file
>> myfile.saveBGF("filename.bgf")
"""

version = "20150317"
this = os.path.basename(sys.argv[0])
if this == "":
	this = "VASP class"


def veclength(v):
	return sqrt(np.dot(v, v))


def angle(a, b):
	angle = arccos(np.dot(a,b) / (veclength(a)*veclength(b)))
	if np.isnan(angle):
		return 0.0
	return rad2deg(angle)



class VASP:
	_latticeParams = [];
	_latticeVectors = [];
	_inv_latticeVectors = [];
	_scaling_factor = 0.0;
	_filename = "";

	_atom_types = []
	_atom_numbers = []
	_selective_dynamics = False;
	_cartesian = False;
	_coords = []
	_filename = ""
	_bgf = None;

	isTriclinic = None;


	def __init__(self, poscar_file=""):
		if poscar_file == "": return;

		self._filename = os.path.abspath(poscar_file)
		print(this + ": " + self._filename + " file loaded.")
		
		vf = open(self._filename)

		# line 1: atom types
		self._atom_types = vf.readline().split()	###

		# line 2: scaling factor
		self._scaling_factor = float(vf.readline())	###
		if self._scaling_factor < 0: nu.die("negative scaling factors are not supported.")

		# line 3-5: pbc
		for i in range(3):
			line = vf.readline()
			h = [self._scaling_factor * float(i) for i in line.split()]
			self._latticeVectors.append(h)		###

		# line 6 or 7: number of atoms OR atom types
		line = vf.readline()
		try:
			self._atom_numbers = [int(i) for i in line.split()]
		except:
			line = vf.readline()
			self._atom_numbers = [int(i) for i in line.split()]

		if len(self._atom_types) != len(self._atom_numbers): nu.die("Number of atom types and atom numbers mismatch.")

		# line 8: coordinate types
		H = np.array(self._latticeVectors, float)
		switch = vf.readline()[0].lower()
		if switch == "s":
			self._selective_dynamics = True
			switch = vf.readline()[0].lower()
		if switch in ("c", "k"):
			self._cartesian = True

		# coordinates
		for index, i in enumerate(self._atom_numbers):
			for j in xrange(i):
				s = vf.readline().split()
				raw_pos = (float(s[0]), float(s[1]), float(s[2]))
				if self._cartesian:
					pos = self._scaling_factor * np.array(raw_pos)
				else:
					pos = np.dot(raw_pos, H)

				self._coords.append(pos)

		self._latticeParams = self.computeLatticeParam(self._latticeVectors)
		print(this + ": Lattice Parameters: " + str(self._latticeParams))

		if self._latticeParams[3] == 90.0 and self._latticeParams[4] == 90.0 and self._latticeParams[5] == 90.0:
			self.isTriclinic = False
		else:
			self.isTriclinic = True
			print(this + ": found triclinic pbc.")

		self._inv_latticeVectors = np.linalg.inv(self._latticeVectors)
		self.convertBGF()


	def computeLatticeParam(self, lv):
		if not len(lv) == 3: nu.die("Proper lattice vectors are not provided.")

		x = lv[0]
		y = lv[1]
		z = lv[2]
		A, B, C = [veclength(v) for v in lv]
		alpha = angle(y, z)
		beta = angle(x, z)
		gamma = angle(x, y)

		return [A, B, C, alpha, beta, gamma]


	def computeLatticeVectors(self, lp):
		# reference: https://pythonhosted.org/MDAnalysis/_modules/MDAnalysis/coordinates/core.html#triclinic_vectors
		B = np.zeros((3,3), dtype=np.float32)
		x, y, z, a, b, c = lp[:6]

		if np.all(lp[:3] == 0):
			return B

		B[0][0] = x
		if a == 90. and b == 90. and c == 90.:
			B[1][1] = y
			B[2][2] = z
		else:
			a = deg2rad(a)
			b = deg2rad(b)
			c = deg2rad(c)
			B[1][0] = y*cos(c)
			B[1][1] = y*sin(c)
			B[2][0] = z*cos(b)
			B[2][1] = z*(cos(a)-cos(b)*cos(c))/sin(c)
			B[2][2] = sqrt(z*z-B[2][0]**2-B[2][1]**2)
		return B


	def distance(self, x1, x2, h, hinv):
		"""
		calculates distance under triclinic pbc
		"""

		s1 = np.dot(hinv, x1); s2 = np.dot(hinv, x2)
		s21 = s2 - s1
		s21 -= np.rint(s21)
		r21 = np.dot(h, s21)
		
		return veclength(r21)


	def convertBGF(self):
		myBGF = bgf.BgfFile()
		natom = 0;

		for index, i in enumerate(self._atom_numbers):
			name = self._atom_types[index]
			name = raw_input("- assign the force field type for " + name + " (default: " + name + ") ?? ") or self._atom_types[index]

			for j in xrange(i):
				atom = bgf.BgfAtom()
				atom.x = self._coords[natom][0]
				atom.y = self._coords[natom][1]
				atom.z = self._coords[natom][2]
				atom.aTag = 1
				atom.aNo = natom + 1
				atom.aName = self._atom_types[index] + str(j + 1)
				atom.rName = "RES"; atom.rNo = 100; atom.chain = "A"
				atom.ffType = name
				myBGF.addAtom(atom)

				natom += 1;

		# assign remarks
		myBGF.BIOGRF = 200
		myBGF.DESCRP = "From " + self._filename + " at " + time.asctime(time.gmtime()) + " on " + os.environ["HOSTNAME"] + " by " + os.environ["USER"]
		myBGF.PERIOD = "111"
		myBGF.AXES = "ZYX"
		myBGF.SGNAME = "P 1                  1    1"
		myBGF.CELLS = [-1, 1, -1, 1, -1, 1]
		myBGF.CRYSTX = self._latticeParams

		self._bgf = myBGF


	def saveBGF(self, *args):
		if self._bgf == None:
			nu.warn("No BGF model converted from " + self._filename + ": save failed.")
			return 0;

		if len(args) != 1:
			tempfile = os.path.abspath('.') + "/" + "POSCAR.bgf"
			nu.warn("Output filename not specified. " + self._filename +" will be stored at " + tempfile)
			self._bgf.saveBGF(tempfile)
		else:
			print(self._filename + " file saved to " + str(args[0]))
			self._bgf.saveBGF(args[0]);


	def find_S_bonds(self):
		if self._bgf == None:
			nu.warn("Cannot make bonds--BGF model unassigned.")
			return 0;

		aNo_pair = [];
		for atom1 in self._bgf.a:
			if not "S" in atom1.ffType: continue;
			x1 = np.array([atom1.x, atom1.y, atom1.z])
			for atom2 in self._bgf.a:
				if atom1 == atom2: continue;
				if not "S" in atom2.ffType: continue;

				x2 = np.array([atom2.x, atom2.y, atom2.z])
				dist = self.distance(x1, x2, self._latticeVectors, self._inv_latticeVectors)
				if 1.95 <= dist <= 2.15:
					aNo_pair.append([atom1.aNo, atom2.aNo])
					
		for i in aNo_pair:
			self._bgf.connectAtoms(i[0], i[1])

		print(this + ": " + str(len(aNo_pair)) + " bonds are created.")
		return 1;


	def make_lammps_data(self, out_file="data.lammps"):
		"""
		Function call:
		>> a.make_lammps_data()	-> creates data.lammps
		>> a.make_lammps_data("data.reax") -> creates data.reax

		This function will create LAMMPS data from VASP instance.
		Requires ReaxFF filename or directory name containing ffield.reax.
		"""
		
		# if no BGF file loaded, LAMMPS data cannot be created.
		if self._bgf == None:
			nu.warn("Cannot make lammps data file--BGF model unassigned.")
			return 0;

		# parse ffield.reax if exists in the specified directory
		reaxfs = "./ffield.reax"
		reaxfs = raw_input("- specify the ReaxFF force field file (default: " + os.path.abspath(reaxfs) + ") ?? ") or reaxfs

		mass = dict()
		if os.path.isdir(reaxfs):
			reaxfs += "/ffield.reax"

		if os.path.isfile(reaxfs):
			reaxff = open(reaxfs)
			while 1:
				line = reaxff.readline()
				if not line: break

				line = line.split()
				for type in self._atom_types:
					if type in line[0]:
						mass[type] = float(line[3])
		else:
			nu.warn("No ffield.reax file found: Abort.")
			return 0

		print(this + ": " + "Sucessfuly loaded ReaxFF file: " + reaxfs)

		# write data.lammps
		fs = open(out_file, 'w')
		fs.write("LAMMPS data file from " + self._filename + " at " + time.asctime(time.gmtime()) + " on " + os.environ["HOSTNAME"] + " by " + os.environ["USER"] + "\n")
		fs.write("\n")
		fs.write("{0:>10}".format(len(self._bgf.a)) + " atoms" + "\n")
		fs.write("\n")
		fs.write("{0:>10}".format(len(self._atom_types)) + " atom types\n")
		fs.write("\n")
		fs.write("\t{0:10.6f} {1:10.6f}".format(0.0, self._latticeParams[0]) + " xlo xhi" + "\n")
		fs.write("\t{0:10.6f} {1:10.6f}".format(0.0, self._latticeParams[1]) + " ylo yhi" + "\n")
		fs.write("\t{0:10.6f} {1:10.6f}".format(0.0, self._latticeParams[2]) + " zlo zhi" + "\n")

		# write xy yz xz if triclinic cell
		if self.isTriclinic:
			l = self._latticeParams
			a = l[0]; b = l[1]; c = l[2]
			cosA = cos(deg2rad(l[3])); cosB = cos(deg2rad(l[4])); cosC = cos(deg2rad(l[5]))
			xy = b * cosC; xz = c* cosB; yz = (b * c * cosA - xy * xz)/b
			fs.write("\t{0:10.6f} {1:10.6f} {2:10.6f} xy xz yz\n".format(xy, xz, yz))
			
		fs.write("\n")
		fs.write("Masses\n")
		fs.write("\n")
		for index, i in enumerate(self._atom_types):
			fs.write("{0:>10}\t{1:10.6f} # {2:<3}\n".format(index + 1, mass[i], i))

		fs.write("\n")
		fs.write("Atoms\n")
		fs.write("\n")
		
		for atom in self._bgf.a:
			fs.write("{0:10} {1:10} {2:10.6f} {3:10.6f} {4:10.6f} {5:10.6f}\n".format(atom.aNo, self._atom_types.index(atom.ffType) + 1, atom.charge, atom.x, atom.y, atom.z))

		print(this + ": " + "Sucessfuly created " + out_file)
		


if __name__ == "__main__":
	print usage
	sys.exit(0)

