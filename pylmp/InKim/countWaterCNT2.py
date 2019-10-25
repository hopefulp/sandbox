#!/home/noische/program/python27/bin/python

import sys, re, string, getopt, optparse, math, time, pprint 
from os import popen
import bgf
import bgftools
import numpy
import operator
import copy
import nutils as nu

import dump

option = ""; args = ""; bgf_file = ""; trj_file = ""; out_file = ""; atom_type = ""; criteria_distance = 0;
usage = """
countWaterCNT.py -b bgf_file -t trj_file
"""

#-----------------
# calculate the distance of specified atoms in every lammps trajectory file
# this is only used during pei analysis by in kim 
# last update: 2010/12/01
#_________________
def countWaterCNT(bgf_file, trj_file, silent=False):

	# init
	timestep = 0; l_timestep = []

	myBGF = bgf.BgfFile(bgf_file)
	myTRJ = open(trj_file)
	myDAT = open(bgf_file[:-4] + ".count.dat", "w")
	myDAT.write("Time\tHeight\tRadius\tNumber\n")

	# extract aNos of CNT in the BGF file
	aNo_CNT = []
	aNo_MtOH_C = []
	for atom in myBGF.a:
		# Carbons in CNT
		if "CNT" in atom.rName:
			aNo_CNT.append(atom.aNo)

		# Carbons in MeOH
		if "MET" in atom.rName and "C" in atom.aName:
			aNo_MtOH_C.append(atom.aNo)

	N_CNT = len(aNo_CNT)	# the number of CNT atoms

	# for every shot in the trajectory file update BGF and manipulate
	myDUMP = dump.dump(trj_file, 0)	# sequential reading

	while 1:

		aNo_MtOH_C_in_CNT = copy.deepcopy(aNo_MtOH_C)

		# initialize for the moment of inertia and the center of mass calculation
		U = 0; Ut = 0; Uv = 0;
		Ixx = 0; Ixy = 0; Ixz = 0;
		Iyx = 0; Iyy = 0; Iyz = 0;
		Izx = 0; Izy = 0; Izz = 0;
		Mx = 0; My = 0; Mz = 0;

		time = myDUMP.next()
		sys.stdout.write('\r' + "Reading timestep.. " + str(time))
		sys.stdout.flush()
		#print("Timestep: " + str(time))
		if time == -1:
			break;

		nu.shutup(); myDUMP.sort(); nu.say();
		l_timestep = myDUMP.time()	# timesteps are 'appended' incrementally
		atoms = myDUMP.viz(len(l_timestep) - 1)[2]	# atom coordinate info: id,type,x,y,z for each atom as 2d array
		box = myDUMP.viz(len(l_timestep) - 1)[1]	# [xlo, ylo, zlo, xhi, yhi, zhi]
		box = [float(i) for i in box]
		boxsize = [box[3]-box[0], box[4]-box[1], box[5]-box[1]]	# simbox size

		myBGF.CRYSTX = [boxsize[0], boxsize[1], boxsize[2], myBGF.CRYSTX[3], myBGF.CRYSTX[4], myBGF.CRYSTX[5]];

		#print("l_atoms " + str(len(atoms)))

		# update the coordinates in the bgf file
		for atom in atoms:
			bgfatom = myBGF.getAtom(atom[0])
			bgfatom.x = float(atom[2])
			bgfatom.y = float(atom[3])
			bgfatom.z = float(atom[4])

		"""
		# com of CNT
		for atom in myBGF.a:
			if "CNT" in atom.rName:
				Mx += atom.x / N_CNT
				My += atom.y / N_CNT
				Mz += atom.z / N_CNT
		"""

		myBGF = bgftools.periodicMoleculeSort(myBGF)

		#print(Mx, My, Mz)
		# transpose for "all atoms in BGF"
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

		# "myBGF" transform
		for atom in myBGF.a:
			v = numpy.matrix([atom.x, atom.y, atom.z]).T
			Uv = Ut * v
			atom.x = float(Uv[0])
			atom.y = float(Uv[1])
			atom.z = float(Uv[2])

		# com of CNT
		Mx = 0; My = 0; Mz = 0;
		for atom in myBGF.a:
			if "CNT" in atom.rName:
				Mx += atom.x / N_CNT
				My += atom.y / N_CNT
				Mz += atom.z / N_CNT
		#print(Mx, My, Mz)

		# transpose for "all atoms in BGF"
		for atom in myBGF.a:
			atom.x -= Mx
			atom.y -= My
			atom.z -= Mz

		""" ### THIS PART (UPDATING LAMMPSTRJ) DOES NOT WORK ###

		myDUMP.aselect.all(time)	# all atoms should be selected to be changed
		myDUMP.tselect.one(time)
		coordx = []; coordy = []; coordz = []
		for atom in myBGF.a:
			coordx.append(atom.x)
			coordy.append(atom.y)
			coordz.append(atom.z)
		myDUMP.setv("xu", coordx)
		myDUMP.setv("yu", coordy)
		myDUMP.setv("zu", coordz)
		"""

		"""	### CHECKING ORIGIN by adding an atom at (0, 0, 0) and is okay. 110914 inkim
		origin = bgf.BgfAtom()
		origin.aTag = 0
		origin.aName = "ORIG"
		origin.ffType = "Ar"
		origin.chain = "A"
		origin.x = 0.0
		origin.y = 0.0
		origin.z = 0.0
		origin.rName = "ORG"
		myBGF.addAtom(origin)
		"""


		# for CNT atoms, calculate some properties
		min_x_CNT = 0.0; max_x_CNT = 0.0; radius_CNT = 0.0; height_CNT = 0.0;

		for aNo in aNo_CNT:
			atom = myBGF.getAtom(aNo)
			# minimum and maximum x coord of CNT: this will be the height of CNT
			if atom.x < min_x_CNT:
				min_x_CNT = atom.x
			if atom.x > max_x_CNT:
				max_x_CNT = atom.x
			radius_CNT += math.sqrt(atom.y**2 + atom.z**2)	# average radius of CNT

		radius_CNT = radius_CNT / N_CNT
		height_CNT = max_x_CNT - min_x_CNT

		# determine whether C in MtOH is in the CNT Cylinder
		for aNo in aNo_MtOH_C:
			atom = myBGF.getAtom(aNo)
			dist = math.sqrt(atom.y**2 + atom.z**2)
			if atom.x > min_x_CNT \
				and atom.x < max_x_CNT \
				and dist < radius_CNT:
				pass;
			else:
				aNo_MtOH_C_in_CNT.remove(aNo)

		#print("Height of CNT: " + str(height_CNT) + " " + str(min_x_CNT) +  " " + str(max_x_CNT))
		#print("Average radius of CNT: " + str(radius_CNT))
		#print("Average number of methanols in CNT: " + str(len(aNo_MtOH_C_in_CNT)))
		myDAT.write(str(time) + "\t")
		myDAT.write(str("{0:<8.3f}".format(height_CNT)))
		myDAT.write(str("{0:<8.3f}".format(radius_CNT)))
		myDAT.write(str(len(aNo_MtOH_C_in_CNT)) + "\n")

		# write output
		myBGF.saveBGF(bgf_file[:-4] + "." + str(time) + ".bgf")

	#myDUMP.write(trj_file.split(".")[0] + "_mod." + trj_file.split(".")[1])
	
	print('')
	return 1

	### end of function


if __name__ == "__main__":
	bgf_file = ""; trj_file = ""

	options, args = getopt.getopt(sys.argv[1:], 'hb:t:a:o:c:', ['help', 'bgf=', 'trj=', 'atom=', 'out=', 'criteria='])

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
	        elif option in ('-t', '--trj'):
	                trj_file = value
	        elif option in ('-a', '--atom'):
	                atom_type = value
	        elif option in ('-o', '--out'):
	                out_file = value
		elif option in ('-c', '--criteria'):
			criteria_distance = float(value)
	        elif option == NULL:
			print(usage)
			sys.exit(0)
	
	# more options: resname CNT, resname solvent, save BGF

	# main call
	countWaterCNT(bgf_file, trj_file)
