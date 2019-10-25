#!/home/noische/program/python27/bin/python

import sys, re, string, getopt, optparse, math, time, pprint 
from os import popen
import bgf
import bgftools
import numpy
import operator

option = ""; args = ""; bgf_file = ""; trj_file = ""; out_file = ""; atom_type = ""; criteria_distance = 0;
usage = """
countWaterCNT.py -b bgf_file -t trj_file -a atom -o out_file -c distance_criteria
"""

# sort keyword
def mysort(x):
	return int(x[0])

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

	# extract aNos in BGF file
	aNo_CNT = []
	for atom in myBGF.a:
		if "CNT" in atom.rName:
			aNo_CNT.append(atom.aNo)

	# for every shot in the trajectory file update BGF and manipulate
	while 1:
		line = myTRJ.readline()
		if not line:
			break

		if 'ITEM: TIMESTEP' in line:
			line = myTRJ.readline()
			timestep = int(line)
			l_timestep.append(timestep)
			sys.stdout.write('\r' + "Reading timestep.." + str(timestep))
			sys.stdout.flush()
			#print("Timestep: " + str(timestep))
			continue

		if 'ITEM: NUMBER OF ATOMS' in line:
			line = myTRJ.readline()
			parse = re.split('\s*', line[:-1])
			n_total_atoms = int(parse[0])
			continue

		if 'ITEM: BOX BOUNDS' in line:
			for i in range(0, 3):
				line = myTRJ.readline()
				line = re.split('\s', line[:-1])
				myBGF.CRYSTX[i] = float(line[1]) - float(line[0])
			continue

		if 'ITEM: ATOMS' in line:
			atoms_location_data = []
			for i in range(0, n_total_atoms):
				line = myTRJ.readline()
				parse = re.split('\s*', line[:-1])
				atom = myBGF.getAtom(int(parse[0]))
				atom.x = float(parse[2]) # update myBGF x coord
				atom.y = float(parse[3]) # update myBGF y coord
				atom.z = float(parse[4]) # update myBGF z coord

		# calculate the moment of inertia (MI) and the center of mass (COM) of CNT from myBGF
		Ixx = 0; Ixy = 0; Ixz = 0;
		Iyx = 0; Iyy = 0; Iyz = 0;
		Izx = 0; Izy = 0; Izz = 0;
		Mx = 0; My = 0; Mz = 0;
		N = len(aNo_CNT)	# the number of CNT atoms

		# com of CNT
		for atom in myBGF.a:
			if "CNT" in atom.rName:
				Mx += atom.x / N
				My += atom.y / N
				Mz += atom.z / N

		# transpose for "all atoms in BGF"
		for atom in myBGF.a:
			atom.x -= Mx
			atom.y -= My
			atom.z -= Mz

		# mi of CNT
		for aNo in aNo_CNT:
			atom = myBGF.getAtom(aNo)
			Ixx += (atom.y**2 + atom.z**2)
			Iyy += (atom.x**2 + atom.z**2)
			Izz += (atom.x**2 + atom.y**2)
			Ixy -= (atom.x * atom.y)
			Ixz -= (atom.x * atom.z)
			Iyz -= (atom.y * atom.z)
			
		I = numpy.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])

		# eigenvalue & eigenvector calculation
		eigval, eigvec = numpy.linalg.eig(I)	# eigval[0] is the minimum among the values.

		# rearrange the U vector
		U = numpy.matrix(eigvec)
		Ut = U.T

		# transform for "all atoms in BGF"
		for atom in myBGF.a:
			v = numpy.matrix([atom.x, atom.y, atom.z]).T
			Uv = Ut * v
			atom.x = Uv[0]
			atom.y = Uv[1]
			atom.z = Uv[2]

		# com of CNT
		Mx = 0; My = 0; Mz = 0;
		for atom in myBGF.a:
			if "CNT" in atom.rName:
				Mx += atom.x / N
				My += atom.y / N
				Mz += atom.z / N

		# transpose for "all atoms in BGF"
		for atom in myBGF.a:
			atom.x -= Mx
			atom.y -= My
			atom.z -= Mz

		# write output
		myBGF.saveBGF(bgf_file[:-4] + "." + str(timestep) + ".bgf")
		### trajectory

	
	print("")
	return 1

	### end of function


if __name__ == "__main__":
	options, args = getopt.getopt(sys.argv[1:], 'hb:t:a:o:c:', ['help', 'bgf=', 'trj=', 'atom=', 'out=', 'criteria='])
	print "Requested options: " + str(options)
	for option, value in options:
	        if option in ('-h', '--help'):
	                print usage
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
			print usage
			sys.exit(0)
	
# main call
countWaterCNT(bgf_file, trj_file)
