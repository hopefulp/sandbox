#!/home/noische/program/python27/bin/python

import sys, re, string, getopt, optparse, math, time, pprint 
from os import popen
import bgf
import bgftools
import numpy
import operator
import copy
import nutils as nu
import itertools
import time

import dump

option = ""; args = ""; bgf_file = ""; trj_file = ""; out_file = ""; 
usage = """
countWaterCNT.py -b bgf_file -t trj_file -f ff_file(s) -n #step -H do_hbond
"""
version = "136193"

#-----------------
# calculate the distance of specified atoms in every lammps trajectory file
# this is only used during pei analysis by in kim 
# last update: 2012/06/01
#_________________
def get_line(trj_file):
	"""
	Open a lazy file stream. This is considered to increase a file reading speed.
	"""
	with open(trj_file) as file:
		for i in file:
			yield i

def sortkey(list):
	return list[0];

def countWaterCNT(bgf_file, trj_file, out_file, silent=False):

	# const
	PI = math.pi
	vdw_r_C = 1.7

	# init
	timestep = 0; l_timestep = []; line = []; n_header = 0;
	t1 = 0; t2 = 0; # clock

	myBGF = bgf.BgfFile(bgf_file)
	myTRJ = open(trj_file)
	myTRJ.seek(0)
	myNewTRJ = open(out_file, 'w')

	# how many steps to go?
	wc_trj_file = popen("grep TIMESTEP " + trj_file + " | wc -l ").read()
	n_timestep = int(wc_trj_file.split()[0]);

	print("The trajectory contains " + str(n_timestep) + " timesteps.")

	# extract aNos of CNT in the BGF file
	aNo_CNT = []
	aNo_MtOH_all = []

	for atom in myBGF.a:
		# Carbons in CNT
		if "CNT" in atom.rName:
			aNo_CNT.append(atom.aNo)

	N_CNT = len(aNo_CNT)	# the number of CNT atoms

	# check if there exists MtOH properly
	if len(aNo_CNT) == 0:
		nu.die("No CNT molecules in the BGF file.")

	# Find header of the trajectory file
	while 1:
		templine = myTRJ.readline()
		line.append(templine.strip('\n').strip('ITEM: '))
		n_header += 1
		if "ITEM: ATOMS" in templine:
			break;

	# INITIAL trajectory information
	timestep = int(line[1])
	natoms = int(line[3])
	boxsize = [line[5].split(' ')[0], line[5].split(' ')[1], line[6].split(' ')[0], line[6].split(' ')[1], line[7].split(' ')[0], line[7].split(' ')[1]]
	boxsize = [float(i) for i in boxsize]
	keywords = line[8].strip('ATOMS ')

	# for every shot in the trajectory file update BGF and manipulate
	dumpatom = get_line(trj_file)
	processed_step = 0;

	t1 = t2 = 0; elapsed_time = 0;

	d_atomtype = dict();	# dict for writing atom types

	# check flag for atom informations: the first shot will be treated, and will be skipped after.
	firsttimeFlag = 0;

	while 1:
		### Show progress
		t1 = time.time();
		remaining_time = elapsed_time * (n_timestep - processed_step)
		sys.stdout.write('\r' + "Reading timestep.. " + str(timestep) + " (Remaining time: " + "{0:4.1f}".format(remaining_time) + " seconds, " + "{0:4.1f} minutes".format(remaining_time/60) + ")")
		sys.stdout.flush()

		processed_step += 1;


		### Read
		try:
			chunk = [next(dumpatom) for i in range(natoms+n_header)]
		except StopIteration:
			break;

		timestep = int(chunk[1])
		natoms = int(chunk[3])
		boxsize = [chunk[5].split(' ')[0], chunk[5].split(' ')[1], chunk[6].split(' ')[0], chunk[6].split(' ')[1], chunk[7].split(' ')[0], chunk[7].split(' ')[1]]; 
		boxsize = [float(i) for i in boxsize]; boxsize = [(boxsize[1] - boxsize[0]), (boxsize[3] - boxsize[2]), (boxsize[5] - boxsize[4])]
		keywords = chunk[8].split('ATOMS ')[1].strip('\n').split(' ')


		### update myBGF with trajectory information ###
		natom_bgf = len(myBGF.a)	# number of atoms in BGF file
	
		if not natom_bgf == natoms:
			nu.die("Number of atoms in trajectory file does not match with BGF file.")
	
		mode = ""
		if 'xs' in keywords or 'ys' in keywords or 'zs' in keywords:
			mode = 'scaled'
		elif 'x' in keywords or 'y' in keywords or 'z' in keywords:
			mode = 'normal'
		elif 'xu' in keywords or 'yu' in keywords or 'zu' in keywords:
			mode = 'unwrapped'
	
		# actual coordinate
		coordinfo = chunk[9:]
	
		# assume that coordinfo is similar to ['id', 'type', 'xs', 'ys', 'zs', 'ix', 'iy', 'iz']
		for atomline in coordinfo:
			atomcoord = atomline.split(' ')
			atom = myBGF.getAtom(int(atomcoord[0]))

			if mode == 'scaled':
				atom.x = float(atomcoord[2]) * boxsize[0]
				atom.y = float(atomcoord[3]) * boxsize[1]
				atom.z = float(atomcoord[4]) * boxsize[2]
	
			elif mode == 'unwrapped':
				atom.x = float(atomcoord[2])
				atom.y = float(atomcoord[3])
				atom.z = float(atomcoord[4])
	
			elif mode == 'normal':
				try:
					ix_index = keywords.index('ix')
					iy_index = keywords.index('iy')
					iz_index = keywords.index('iz')
				except ValueError:
					nu.warn("No image information no the trajectory file. Will be treated as unwrapped.")
					atom.x = float(atomcoord[2])
					atom.y = float(atomcoord[3])
					atom.z = float(atomcoord[4])
				else:
					atom.x = (int(atomcoord[ix_index]) * boxsize[0]) + float(atomcoord[2]) 
					atom.y = (int(atomcoord[iy_index]) * boxsize[1]) + float(atomcoord[3]) 
					atom.z = (int(atomcoord[iz_index]) * boxsize[2]) + float(atomcoord[4]) 
	
			try:
				for i in range(0, 3):
					myBGF.CRYSTX[i] = boxsize[i]
			except:
				pass;
				#nu.warn("Crystal information error: is this file not periodic?")
				
			### for the first shot, record atom types
			if firsttimeFlag == 0:
				d_atomtype[int(atomcoord[0])] = atomcoord[1]
		
		### myBGF update complete! ###

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
			atom.x += Mx
			atom.y += My
			atom.z += Mz

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

		### update crystx information: note that CNT is along x-axis.	# is this right?? 20130619
		if boxsize[0] > boxsize[1] and boxsize[0] > boxsize[2]:
			boxsize = [boxsize[0], boxsize[1], boxsize[2]]
		elif boxsize[1] > boxsize[0] and boxsize[1] > boxsize[2]:
			boxsize = [boxsize[1], boxsize[0], boxsize[2]]
		elif boxsize[2] > boxsize[0] and boxsize[2] > boxsize[1]:
			boxsize = [boxsize[2], boxsize[0], boxsize[1]]

		for i in range(0, 3):
			myBGF.CRYSTX[i] = boxsize[i]

		# transpose to the center of box
		for atom in myBGF.a:
			atom.x += (boxsize[0] / 2)
			atom.y += (boxsize[1] / 2)
			atom.z += (boxsize[2] / 2)

		# apply periodic condition
		myBGF = bgftools.periodicMoleculeSort(myBGF, 0)

		### Write LAMMPS trjectory
		output = "ITEM: TIMESTEP\n"
		output += str(timestep) + "\n"
		output += "ITEM: NUMBER OF ATOMS\n"
		output += str(natom_bgf) + "\n"
		output += "ITEM: BOX BOUNDS pp pp pp\n"
		output += "0 " + str(boxsize[0]) + "\n"
		output += "0 " + str(boxsize[1]) + "\n"
		output += "0 " + str(boxsize[2]) + "\n"
		output += "ITEM: ATOMS id type xu yu zu\n"
		for atom in myBGF.a:
			output += str(atom.aNo) + " " + d_atomtype[atom.aNo] + " " + str(atom.x) + " " + str(atom.y) + " " + str(atom.z) + "\n"

		myNewTRJ.write(output)
	
		t2 = time.time()	# time mark
		elapsed_time = t2 - t1;

	myNewTRJ.close()
	print('')

	return 1

	### end of function


if __name__ == "__main__":
	bgf_file = ""; trj_file = ""; out_file = "";

	options, args = getopt.getopt(sys.argv[1:], 'hb:t:o:', ['help', 'bgf=', 'trj=', 'out='])

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
		elif option in ('-o', '--out'):
			out_file = value
	        elif option == NULL:
			print(usage)
			sys.exit(0)
	
	# more options: resname CNT, resname solvent, save BGF

	# main call
	countWaterCNT(bgf_file, trj_file, out_file)
