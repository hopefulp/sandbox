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
import cPickle as pkl

option = ""; args = ""; bgf_file = ""; trj_file = ""; out_file = ""; 
usage = """
countWaterCNT.py -b bgf_file -t trj_file -n #step -s avg_step
"""
version = "130726"

interval = 0.1;

#-----------------
# This calculates force on water molecules: method 2
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

def getSlipLength(bgf_file, trj_file, ff_file, n_step, n_avg_step, silent=False):

	# const
	PI = math.pi
	vdw_r_C = 1.7

	# init
	timestep = 0; l_timestep = []; line = []; n_header = 0;
	t1 = 0; t2 = 0; # clock

	l_data = [];	# stores [r vz]
	l_F_N = []; l_F_N_inner = []; l_F_N_outer = [];	# stores sum of friction force
	l_result = [];
	l_avg_radius_CNT = [];

	n_avg_count = 0;
	n_avg_count_step = 0;

	myBGF = bgf.BgfFile(bgf_file)
	myTRJ = open(trj_file)
	myTRJ.seek(0)
	myOUT = open(trj_file + ".fz.profile", 'w')


	# how many steps to go?
	if not silent: print("== Step 1. Loading LAMMPS Trajectory")
	wc_trj_file = popen("grep TIMESTEP " + trj_file + " | wc -l ").read()
	n_timestep = int(wc_trj_file.split()[0]);
	if n_step == 0:
		n_step = n_timestep;

	print("The trajectory contains " + str(n_timestep) + " timesteps.")

	# extract aNos of CNT in the BGF file
	aNo_CNT = []
	aNo_WAT_O = []

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
	keywords = line[8].strip('ATOMS ')	# assume ATOMS id type xu yu zu vx vy vz fx fy fz

	# for every shot in the trajectory file update BGF and manipulate
	dumpatom = get_line(trj_file)
	processed_step = 0;

	t1 = t2 = 0; elapsed_time = 0;

	while 1:
		### Show progress
		#t1 = time.time();
		#remaining_time = elapsed_time * (n_step - processed_step)
		sys.stdout.write('\r' + "Reading timestep.. " + str(timestep))# + " (Remaining time: " + "{0:4.1f}".format(remaining_time) + " seconds, " + "{0:4.1f} minutes".format(remaining_time/60) + ")")
		sys.stdout.flush()

		if processed_step == n_step:
			break;

		### Read
		try:
			chunk = [next(dumpatom) for i in range(natoms+n_header)]
		except StopIteration:
			break;

		timestep = int(chunk[1])
		l_timestep.append(timestep)
		natoms = int(chunk[3])
		boxsize = [chunk[5].split(' ')[0], chunk[5].split(' ')[1], chunk[6].split(' ')[0], chunk[6].split(' ')[1], chunk[7].split(' ')[0], chunk[7].split(' ')[1]]; 
		boxsize = [float(i) for i in boxsize]; boxsize = [(boxsize[1] - boxsize[0]), (boxsize[3] - boxsize[2]), (boxsize[5] - boxsize[4])]
		keywords = chunk[8].split('ATOMS ')[1].strip('\n').split(' ')


		### update myBGF with trajectory information ###
		natom_bgf = len(myBGF.a)	# number of atoms in BGF file
	
		if not natom_bgf == natoms:
			nu.warn("Number of atoms in trajectory file does not match with BGF file.")
	
		mode = 'unwrapped'	# assume that "dump            1 all custom 100 ${sname}${rtemp}K.nvt.lammps id type xu yu zu vx vy vz fx fy fz" in lammps input
	
		# actual coordinate
		coordinfo = chunk[9:]
	
		# modified for fast treatment
		for atomline in coordinfo:
			atomcoord = atomline.split(' ')
			atom = myBGF.getAtom(int(atomcoord[0]))
	
			atom.x = float(atomcoord[2])	# coord
			atom.y = float(atomcoord[3])
			atom.z = float(atomcoord[4])
			atom.vx = float(atomcoord[5])	# vel
			atom.vy = float(atomcoord[6])
			atom.vz = float(atomcoord[7])
			if 'fx' in keywords:
				atom.fx = float(atomcoord[8])	# force
				atom.fy = float(atomcoord[9])
				atom.fz = float(atomcoord[10])

			try:
				for i in range(0, 3):
					myBGF.CRYSTX[i] = boxsize[i]
			except:
				pass;
				#nu.warn("Crystal information error: is this file not periodic?")
				
		# apply periodic condition
		#myBGF = bgftools.periodicMoleculeSort(myBGF, 0, True)

		### myBGF update complete! ###


		# calculate sum of tangential forces on carbon atoms of CNT (sum fz)
		F_N = 0.0; 
		for aNo in aNo_CNT:
			atom = myBGF.getAtom(aNo)
			F_N += atom.fz

		l_F_N.append([timestep, F_N]);


		# ending for timestep treatment
		#t2 = time.time()	# time mark
		#elapsed_time = t2 - t1;
		processed_step += 1


	output = "";
	for i in l_F_N:
		output += str(i[0]) + '\t' + str(i[1]) + '\n'

	myOUT.write(output)


	print('')
	return 1

	### end of function


if __name__ == "__main__":
	# n_avg_step = how many steps will be averaged
	bgf_file = ""; trj_file = ""; ff_file = ""; n_step = 0; n_avg_step = 1; ff_file = ""; 

	options, args = getopt.getopt(sys.argv[1:], 'hb:t:n:s:f:', ['help', 'bgf=', 'trj=', 'number=', 'step=', 'ff='])

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
		elif option in ('-n', '--number'):
			n_step = int(value)
		elif option in ('-s', '--step'):
			n_avg_step = int(value)
		elif option in ('-f', '--ff'):
			ff_file = value
	        elif option == NULL:
			print(usage)
			sys.exit(0)
	
	# more options: resname CNT, resname solvent, save BGF

	# main call
	getSlipLength(bgf_file, trj_file, ff_file, n_step, n_avg_step)
