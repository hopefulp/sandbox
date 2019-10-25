#!/home/noische/program/python27/bin/python

import sys, re, string, getopt, optparse, math, time, pprint 
import bgf
import bgftools
import numpy
import operator
import copy
import nutils as nu

option = ""; args = ""; bgf_file = ""; 
usage = """
CNT_alignCNT.py -b bgf_file -f ff_file(s)
"""
version = "130615"

"""
20130615 there is a problem that crystx information does not be updated
"""
def alignCNT(bgf_file, ff_file, silent=False):

	# const
	PI = math.pi
	vdw_r_C = 1.7

	# init

	myBGF = bgf.BgfFile(bgf_file)

	# extract aNos of CNT in the BGF file
	aNo_CNT = []
	aNo_MtOH_C = []
	aNo_MtOH_all = []

	for atom in myBGF.a:
		# Carbons in CNT
		if "CNT" in atom.rName:
			aNo_CNT.append(atom.aNo)

		# Carbons in MeOH
		if "MET" in atom.rName and "C" in atom.aName:
			aNo_MtOH_C.append(atom.aNo)

	N_CNT = len(aNo_CNT)	# the number of CNT atoms

	# check if there exists MtOH properly
	#if len(aNo_MtOH_C) == 0:
	#	nu.die("No MtOH molecules in the BGF file.")
	if len(aNo_CNT) == 0:
		nu.die("No CNT molecules in the BGF file.")


	# from here @@#$!@$!T$
	aNo_MtOH_C_atoms = []
	aNo_atoms_in_CNT = []

	aNo_MtOH_C_in_CNT = copy.deepcopy(aNo_MtOH_C)

	### Realign whole system to x axis
	sys.stdout.write(' Realigning.. ')
	sys.stdout.flush()
	# initialize for the moment of inertia and the center of mass calculation
	U = 0; Ut = 0; Uv = 0;
	Ixx = 0; Ixy = 0; Ixz = 0;
	Iyx = 0; Iyy = 0; Iyz = 0;
	Izx = 0; Izy = 0; Izz = 0;
	Mx = 0; My = 0; Mz = 0;

	# transpose for "all atoms in BGF": move COM of CNT as origin
	# com of CNT
	Mx = 0; My = 0; Mz = 0;
	for atom in myBGF.a:
		if "CNT" in atom.rName:
			Mx += atom.x / N_CNT
			My += atom.y / N_CNT
			Mz += atom.z / N_CNT

	# move
	for atom in myBGF.a:
		atom.x -= Mx
		atom.y -= My
		atom.z -= Mz

	# calculate the moment of inertia (MI) and the center of mass (COM) of CNT from "myBGF"
	for aNo in aNo_CNT:
		atom = myBGF.getAtom(aNo)
		Ixx += (atom.y**2 + atom.z**2) / N_CNT
		Iyy += (atom.x**2 + atom.z**2) / N_CNT
		Izz += (atom.x**2 + atom.y**2) / N_CNT
		Ixy -= (atom.x * atom.y) / N_CNT
		Ixz -= (atom.x * atom.z) / N_CNT
		Iyz -= (atom.y * atom.z) / N_CNT
		
	# the moment of inertia tensor
	I = numpy.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])

	# eigenvalue & eigenvector calculation
	eigval, eigvec = numpy.linalg.eig(I)	# eigval[0] is the minimum among the values.

	# rearrange the U vector
	U = numpy.matrix(eigvec)
	Ut = U.T

	# "myBGF" rotation
	for atom in myBGF.a:
		v = numpy.matrix([atom.x, atom.y, atom.z]).T
		Uv = Ut * v
		atom.x = float(Uv[0])
		atom.y = float(Uv[1])
		atom.z = float(Uv[2])

	# for CNT atoms, calculate some properties
	min_x_CNT = 0.0; max_x_CNT = 0.0; radius_CNT = 0.0; height_CNT = 0.0; aNo_MtOH_C_not_in_CNT = []; temp = [];
	min_y_CNT = 0.0; max_y_CNT = 0.0;
	min_z_CNT = 0.0; max_z_CNT = 0.0;
	CNT_orientation = "";

	# check the orientation of CNT
	for aNo in aNo_CNT:
		atom = myBGF.getAtom(aNo)
		# minimum and maximum x coord of CNT: this will be the height of CNT
		if atom.x < min_x_CNT:
			min_x_CNT = atom.x
		if atom.x > max_x_CNT:
			max_x_CNT = atom.x
		if atom.y < min_y_CNT:
			min_y_CNT = atom.y
		if atom.y > max_y_CNT:
			max_y_CNT = atom.y
		if atom.z < min_z_CNT:
			min_z_CNT = atom.z
		if atom.z > max_z_CNT:
			max_z_CNT = atom.z

	x_diff = max_x_CNT - min_x_CNT
	y_diff = max_y_CNT - min_y_CNT
	z_diff = max_z_CNT - min_z_CNT

	if x_diff > y_diff and x_diff > z_diff:
		# CNT is aligned along x axis
		height_CNT = x_diff
		CNT_orientation = "x"
	elif y_diff > x_diff and y_diff > z_diff:
		# CNT is aligned along y axis
		height_CNT = y_diff
		CNT_orientation = "y"
	elif z_diff > x_diff and z_diff > y_diff:
		# CNT is aligned along z axis
		height_CNT = z_diff
		CNT_orientation = "z"

	if CNT_orientation == "x":
		for aNo in aNo_CNT:
			atom = myBGF.getAtom(aNo)
			radius_CNT += math.sqrt(atom.y**2 + atom.z**2)	# average radius of CNT
	elif CNT_orientation == "y":
		for aNo in aNo_CNT:
			atom = myBGF.getAtom(aNo)
			radius_CNT += math.sqrt(atom.x**2 + atom.z**2)	# average radius of CNT
	elif CNT_orientation == "z":
		for aNo in aNo_CNT:
			atom = myBGF.getAtom(aNo)
			radius_CNT += math.sqrt(atom.x**2 + atom.y**2)	# average radius of CNT

	radius_CNT = radius_CNT / N_CNT

	### Rotate y or z to x axis
	Tyx = numpy.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
	Tzx = numpy.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])

	if CNT_orientation == "y":
		for atom in myBGF.a:
			u = numpy.matrix([atom.x, atom.y, atom.z]).T
			Tu = Tyx*u
			atom.x = float(Tu[0])
			atom.y = float(Tu[1])
			atom.z = float(Tu[2])
			min_x_CNT = min_y_CNT;
			max_x_CNT = max_y_CNT;

	elif CNT_orientation == "z":
		for atom in myBGF.a:
			u = numpy.matrix([atom.x, atom.y, atom.z]).T
			Tu = Tzx*u
			atom.x = float(Tu[0])
			atom.y = float(Tu[1])
			atom.z = float(Tu[2])
			min_x_CNT = min_z_CNT;
			max_x_CNT = max_z_CNT;

	# write output
	myBGF.saveBGF(bgf_file[:-4] + ".aligned.bgf")

	return 1

	### end of function


if __name__ == "__main__":
	bgf_file = ""; ff_file = ""; 
	options, args = getopt.getopt(sys.argv[1:], 'hb:f:', ['help', 'bgf=', 'forcefield='])

	if len(sys.argv) < 2:
		print(usage)
		sys.exit(0)

	print "Requested options: " + str(options)

	for option, value in options:
	        if option in ('-h', '--help'):
	                print(usage)
			sys.exit(0)
	        elif option in ('-b', '--bgf'):
	                bgf_file = value
	        elif option in ('-f', '--forcefield'):
	                ff_file = str(value).strip()
	        elif option == NULL:
			print(usage)
			sys.exit(0)
	
	# main call
	alignCNT(bgf_file, ff_file)
