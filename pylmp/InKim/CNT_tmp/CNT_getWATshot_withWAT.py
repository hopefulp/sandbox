#!/home/noische/program/python27/bin/python

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

option = ""; args = ""; bgf_file = ""; trj_file = ""; out_file = ""; 
usage = """
countWaterCNT.py -b bgf_file -t trj_file -n #step
"""
version = "130717"

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
	interval = 0.1;

	l_data = [];	# stores vz
	l_avg_radius_CNT = [];

	myBGF = bgf.BgfFile(bgf_file)
	myTRJ = open(trj_file)
	myTRJ.seek(0)

	# how many steps to go?
	wc_trj_file = popen("grep TIMESTEP " + trj_file + " | wc -l ").read()
	n_timestep = int(wc_trj_file.split()[0]);

	print("The trajectory contains " + str(n_timestep) + " timesteps.")

	# extract aNos of CNT in the BGF file
	aNo_CNT = []
	aNo_WAT_O = []
	aNo_WAT_all = []

	for atom in myBGF.a:
		# Carbons in CNT
		if "CNT" in atom.rName:
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
		sys.stdout.write('\r' + "Reading timestep.. " + str(timestep) )
		sys.stdout.flush()

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

		if timestep != n_step:
			continue;

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

			try:
				for i in range(0, 3):
					myBGF.CRYSTX[i] = boxsize[i]
			except:
				pass;
				#nu.warn("Crystal information error: is this file not periodic?")
				
		# apply periodic condition
		myBGF = bgftools.periodicMoleculeSort(myBGF, 0, True)

		### myBGF update complete! ###

		#aNo_WAT_O_in_CNT = copy.deepcopy(aNo_WAT_O)


		### align CNT to z axis
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
			atom.x = float(Uv[1])
			atom.y = float(Uv[2])
			atom.z = float(Uv[0])


		# for CNT atoms, calculate some properties
		min_x_CNT = 1000.0; max_x_CNT = -1000.0; radius_CNT = 0.0; height_CNT = 0.0; temp = [];
		min_y_CNT = 1000.0; max_y_CNT = -1000.0;
		min_z_CNT = 1000.0; max_z_CNT = -1000.0;
		CNT_orientation = "";

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

		# radius of CNT
		for aNo in aNo_CNT:
			atom = myBGF.getAtom(aNo)
			l_radius_CNT.append( math.sqrt( (atom.x - x_CNT)**2 + (atom.y - y_CNT)**2 ))

		radius_CNT, a = nu.meanstdv(l_radius_CNT)
		l_avg_radius_CNT.append(radius_CNT)


		### get water molecules in CNT and out of CNT

		# inside the CNT := min_z_CNT <= z <= max_z_CNT and (x - (x_diff/2))**2 + (y - (y_diff/2))**2 < radius_CNT**2
		# aNo_WAT_O_atoms: molecules which O atom is within CNT
		# we don't need to calculate H atoms. Let's consider only O atoms
		margin = 0.0;	# water molecules far from the margin will be only considered
		aNo_WAT_O_in_CNT = []
		aNo_WAT_O_not_in_CNT = []
		myBGF2 = copy.deepcopy(myBGF)

		# delete graphenes
		for atom in myBGF2.a:
			if "GRA" in atom.rName:
				aNo_WAT_O_not_in_CNT.append(myBGF.a2i[atom.aNo])

		# delete water molecules without the length of CNT
		n_wat_in_CNT = 0;
		n_wat_out_CNT = 0;
		for aNo in aNo_WAT_O:
			atom = myBGF.getAtom(aNo)
			dist_sq = (atom.x - x_CNT)**2 + (atom.y - y_CNT)**2
			# water in CNT
			if "WAT" in atom.rName and "O" in atom.ffType and min_z_CNT + margin <= atom.z and atom.z <= max_z_CNT - margin and dist_sq < radius_CNT**2:
				aNo_WAT_O_in_CNT.append(aNo)
				n_wat_in_CNT += 1;
			# water not in CNT but inside the box
			if "WAT" in atom.rName and atom.z < max_z_CNT and atom.z > min_z_CNT:
				aNo_WAT_O_in_CNT.append(aNo)
				n_wat_out_CNT += 1;
			else:
				delete_list = [];
				dummy = bgftools.getmolecule(myBGF, myBGF.getAtom(aNo), delete_list)
				for i in delete_list:
					aNo_WAT_O_not_in_CNT.append(myBGF.a2i[i])

		myBGF2.delAtoms(aNo_WAT_O_not_in_CNT)
		myBGF2.renumber()
		myBGF2.saveBGF(bgf_file[:-4] + "." + str(timestep) + "." + str(len(aNo_WAT_O_in_CNT)) + ".bgf")

		sys.stdout.write('Done                                                        ')
		sys.stdout.flush()

		#print(len(aNo_WAT_O_in_CNT))
		print(n_wat_in_CNT)
		print(n_wat_out_CNT)

		processed_step += 1;


	print('')

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
