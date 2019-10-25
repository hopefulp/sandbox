#!/home/noische/python

import sys, re, string, getopt, optparse, math, time, pprint 
from os import popen
import bgf
import bgftools
import operator
import copy
import nutils as nu
import itertools
import time
import numpy as np
import lammpstools as lt
import os

option = ""; args = ""; bgf_file = ""; trj_file = ""; out_file = ""; 
usage = """
Count the number of molecules in CNT with LAMMPS trajectory file.

countWaterCNT.py -b bgf_file -t trj_file -n #step
"""
version = "140304"

#-----------------
# Get averaged velocity of water molecules confined in CNT
# 
# 
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

def getVelocity(bgf_file, trj_file, n_step, silent=False):

	# const
	PI = math.pi
	vdw_r_C = 1.7

	# init
	timestep = 0; l_timestep = []; line = []; n_header = 0;
	t1 = 0; t2 = 0; # clock

	l_data = [];	# stores vz
	l_avg_radius_CNT = [];

	myBGF = bgf.BgfFile(bgf_file)
	myTRJ = open(trj_file)
	myTRJ.seek(0)

	ftemp = open("countWater42.profile", 'w')
	ftemp.write(str(sys.argv) + "\n")

        curr_dir = os.path.abspath(".")
        temp_dir = curr_dir + "/CountWAT42/"
        if not os.path.isdir(temp_dir): os.makedirs(temp_dir)

	# how many steps to go?
	n_timestep = len(lt.getTrjInfo(trj_file))
	if n_step == 0:
		n_step = n_timestep;

	print(" ..The trajectory contains " + str(n_timestep) + " timesteps.")
	print("The script will proceed for the last " + str(n_step) + " timesteps.")

	# extract aNos of CNT in the BGF file
	aNo_CNT = []
	aNo_WAT_O = []
	aNo_WAT_all = []

	for atom in myBGF.a:
		# Carbons in CNT or atoms in BNNT
		if "NT" in atom.rName:
			aNo_CNT.append(atom.aNo)

		# Oxygen in water
		if "WAT" in atom.rName and "O" in atom.aName:
			aNo_WAT_O.append(atom.aNo)

	N_CNT = len(aNo_CNT)	# the number of CNT atoms

	# check if there exists water properly
	if len(aNo_WAT_O) == 0:
		nu.die("No water molecules in the BGF file.")
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

	while 1:

		### Show progress
		t1 = time.time();
		remaining_time = elapsed_time * (n_step - processed_step)
		sys.stdout.write('\r' + "Reading timestep.. " + str(timestep) + " (Remaining time: " + "{0:4.1f}".format(remaining_time) + " seconds, " + "{0:4.1f} minutes".format(remaining_time/60) + ")")
		sys.stdout.flush()

		if processed_step == n_step:
			break;

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
	
		mode = 'unwrapped'	# assume that "dump            1 all custom 100 ${sname}${rtemp}K.nvt.lammps id type xu yu zu vx vy vz" in lammps input
	
		# actual coordinate
		coordinfo = chunk[9:]
	
		# modified for fast treatment
		for atomline in coordinfo:
			atomcoord = atomline.split(' ')
			atom = myBGF.getAtom(int(atomcoord[0]))
	
			atom.x = float(atomcoord[2])
			atom.y = float(atomcoord[3])
			atom.z = float(atomcoord[4])
			atom.vx = float(atomcoord[5])
			atom.vy = float(atomcoord[6])
			atom.vz = float(atomcoord[7])

		if myBGF.CRYSTX != []:
			for i in range(0, 3):
				myBGF.CRYSTX[i] = boxsize[i]
		else:
			nu.warn("Crystal information error: is this file not periodic?")
			for i in range(0, 3):
				myBGF.CRYSTX.append(boxsize[i])

		### myBGF update complete! ###

		### align CNT to z axis
		# initialize for the moment of inertia and the center of mass calculation
		U = 0; Ut = 0; Uv = 0;
		Ixx = 0; Ixy = 0; Ixz = 0;
		Iyx = 0; Iyy = 0; Iyz = 0;
		Izx = 0; Izy = 0; Izz = 0;
		Mx = 0; My = 0; Mz = 0;

		# set COM of CNT as origin
		Mx = 0; My = 0; Mz = 0;	# Center of mass of CNT
		for atom in myBGF.a:
			if "NT" in atom.rName:
				Mx += atom.x / N_CNT
				My += atom.y / N_CNT
				Mz += atom.z / N_CNT

		for atom in myBGF.a:
			atom.x -= Mx	# move
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
		I = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])

		# eigenvalue & eigenvector calculation
		eigval, eigvec = np.linalg.eig(I)	# eigval[0] is the minimum among the values.
	
		# rearrange the U vector
		U = np.matrix(eigvec)
		Ut = U.T

		# "myBGF" rotation
		for atom in myBGF.a:
			v = np.matrix([atom.x, atom.y, atom.z]).T
			Uv = Ut * v
			atom.x = float(Uv[2])
			atom.y = float(Uv[1])
			atom.z = float(Uv[0])
		#dimension = np.matrix(boxsize).T
		#boxsize_prime = np.array(np.dot(Ut, dimension).T)
		boxsize_prime = np.transpose(np.dot(Ut,np.transpose(np.array(boxsize))))	# new pbc from jackjack5
		print(boxsize)
		print(boxsize_prime)

		# move atoms to box center
		for atom in myBGF.a:
			atom.x -= boxsize_prime[0]/2
			atom.y -= boxsize_prime[1]/2
			atom.z -= boxsize_prime[2]/2

		myBGF.saveBGF(temp_dir + bgf_file.split(".bgf")[0] + "." + str(timestep) + ".beforepbc.bgf")
		myBGF = bgftools.periodicMoleculeSort(myBGF, 0, myBGF.CRYSTX)
		myBGF.saveBGF(temp_dir + bgf_file.split(".bgf")[0] + "." + str(timestep) + ".bgf")

		# for CNT atoms, calculate some properties
		min_x_CNT = 1000.0; max_x_CNT = -1000.0; radius_CNT = 0.0; height_CNT = 0.0; 
		min_y_CNT = 1000.0; max_y_CNT = -1000.0;
		min_z_CNT = 1000.0; max_z_CNT = -1000.0;

		l_radius_CNT = [];
		l_x_CNT = [];
		l_y_CNT = [];
		for aNo in aNo_CNT:
			atom = myBGF.getAtom(aNo)

			# center of CNT
			l_x_CNT.append(atom.x)
			l_y_CNT.append(atom.y)

			# height
			if atom.z < min_z_CNT:
				min_z_CNT = atom.z
			if atom.z > max_z_CNT:
				max_z_CNT = atom.z
	
		x_CNT, a = nu.meanstdv(l_x_CNT)
		y_CNT, a = nu.meanstdv(l_y_CNT)
		z_diff = max_z_CNT - min_z_CNT
		height_CNT = z_diff

		# radius of CNT
		for aNo in aNo_CNT:
			atom = myBGF.getAtom(aNo)
			l_radius_CNT.append( math.sqrt( (atom.x - x_CNT)**2 + (atom.y - y_CNT)**2 ))

		radius_CNT, a = nu.meanstdv(l_radius_CNT)
		l_avg_radius_CNT.append(radius_CNT)


		### get water molecules in CNT

		# inside the CNT := min_z_CNT <= z <= max_z_CNT and (x - (x_diff/2))**2 + (y - (y_diff/2))**2 < radius_CNT**2
		# aNo_WAT_O_atoms: molecules which O atom is within CNT
		# we don't need to calculate H atoms. Let's consider only O atoms
		margin = 0.0;	# water molecules far from the margin will be only considered
		aNo_WAT_O_in_CNT = []
		for aNo in aNo_WAT_O:
			atom = myBGF.getAtom(aNo)
			dist_sq = (atom.x - x_CNT)**2 + (atom.y - y_CNT)**2
			if "WAT" in atom.rName and "O" in atom.ffType and min_z_CNT + margin <= atom.z and atom.z <= max_z_CNT - margin and dist_sq < radius_CNT**2:
				aNo_WAT_O_in_CNT.append(aNo)
			else:
				pass;

		### record number of water molecules in CNT
		#l_data.append([timestep, len(aNo_WAT_O_in_CNT), radius_CNT, z_diff])

		#output = str(timestep) + '\t' + str(len(aNo_WAT_O_in_CNT)) + '\t' + str(radius_CNT) + '\t' + str(z_diff) + '\n'
		output = str(timestep) + '\t' + str(len(aNo_WAT_O_in_CNT)) + '\t' + str(radius_CNT) + '\t' + str(min_z_CNT) + '\t' + str(max_z_CNT) + '\t' + str(z_diff) + '\t' + str(height_CNT) + '\t' + str(len(aNo_CNT)) + '\n'	# debug
		ftemp.write(output)


		#sys.stdout.write('Done                                                        ')
		sys.stdout.flush()

		t2 = time.time()	# time mark
		elapsed_time = t2 - t1;
		processed_step += 1;


	print('')
	ftemp.close()
	print("numbers of water molecules are written in countWater.profile ..Done.")


	return 1

	### end of function


if __name__ == "__main__":
	bgf_file = ""; trj_file = ""; ff_file = ""; n_step = 0; 

	options, args = getopt.getopt(sys.argv[1:], 'hb:t:n:', ['help', 'bgf=', 'trj=', 'step='])

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
		elif option in ('-n', '--step'):
			n_step = int(value)
	        elif option == NULL:
			print(usage)
			sys.exit(0)
	
	# more options: resname CNT, resname solvent, save BGF

	# main call
	getVelocity(bgf_file, trj_file, n_step)
