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
import lammpstools as lt

option = ""; args = ""; 
A = 0.0; v = 0.0;

usage = """
countWaterCNT.py -b bgf_file -t trj_file -n #step
"""
version = "130729"

#-----------------
# 
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

def getSlipLength(bgf_file, trj_file, force_profile, silent=False):
	"""
	This calculates the profile of averaged velocity according to r. 
	NOTE: This is different from averaged velocity profile according to r.
	"""

	# const
	PI = math.pi
	vdw_r_C = 1.7

	# init
	timestep = 0; l_timestep = []; line = []; n_header = 0; output = "";
	t1 = 0; t2 = 0; # clock
	output += "t\tforce\tarea\tvz\tvisc\tsliplen\n"

	l_radial_vz_profile = [];
	l_avg_data = [];
	l_result = [];

	l_avg_count = [];
	n_avg_count_step = 0;

	myBGF = bgf.BgfFile(bgf_file)
	myTRJ = open(trj_file)
	myTRJ.seek(0)
	myOUT = open(trj_file + ".sliplength.profile", 'w')


	# how many steps to go?
	if not silent: print("Scanning trajectory file..")
	l_timestep = lt.getTrjInfo(trj_file, silent=False)
	n_timestep = len(l_timestep)

	if not silent: print("\nThe trajectory contains " + str(n_timestep) + " timesteps.")


	# read force data from force profile
	if not silent: print("Reading force profile..")
	f_fp = open(force_profile)
	force = dict();
	while 1:
		line = f_fp.readline()
		if not line:
			break;
		if line[0] == "#":
			continue;
		parse = line.split()
		step = int(parse[0])
		force[step] = float(parse[1])


	# extract aNos of CNT in the BGF file
	aNo_CNT = []; aNo_WAT_O = []; aNo_WAT_all = []

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

	t1 = t2 = 0; elapsed_time = 0;

	# Find header of the trajectory file
	line = []; n_header = 0;
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

	# calculate properties
	dumpatom = get_line(trj_file)
	processed_step = 0;
	while 1:
		### Show progress
		t1 = time.time();
		remaining_time = elapsed_time * (n_timestep - processed_step)

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
		sys.stdout.write('\r' + "Fetched timestep: " + str(timestep) + " (Remaining time: " + "{0:4.1f}".format(remaining_time) + " seconds, " + "{0:4.1f} minutes".format(remaining_time/60) + ")")
		sys.stdout.flush()

		if timestep == 0:
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
	
			atom.vx = float(atomcoord[5])
			atom.vy = float(atomcoord[6])
			atom.vz = float(atomcoord[7])

		### myBGF update complete! ###


		# calculate z-directional velocity of water molecules wrt r		vcm = sum vi mi / sum mi
		mass_O = 16.0000
		mass_H = 1.008
		mass_H2O = mass_O + 2 * mass_H	# H2O mass
		vz = [];	# data container

		# z velocity calculation
		for ano in aNo_WAT_O:
			atom = myBGF.getAtom(ano)	# oxygen
			vcmz = atom.vz * mass_O

			for ano2 in atom.CONECT:
				atom2 = myBGF.getAtom(ano2)	# hydrogens
				vcmz += atom2.vz * mass_H

			vcmz /= mass_H2O

			vz.append(vcmz*(10**5))	# velocity in A/fs -> m/s
			
		avg_vz, _ = nu.meanstdv(vz)

		f = force[timestep] * (6.95e9)
		l = f / A / avg_vz	# lambda
		b = v / l * 1e9	# slip length

		output += str(timestep) + "\t" + str(force[timestep]) + "\t" + str(A) + "\t" + str(avg_vz) + "\t" + str(v) + "\t" +  str(b) + "\n"

		t2 = time.time()	# time mark
		elapsed_time = t2 - t1;
		processed_step += 1

	myOUT.write(output)
	myOUT.close()

	print('')
	return 1

	### end of function


if __name__ == "__main__":
	bgf_file = ""; trj_file = ""; n_pass = 0; force_profile = ""; 
	options, args = getopt.getopt(sys.argv[1:], 'hb:t:i:p:A:v:', ['help', 'bgf=', 'trj=', 'force=', 'area=', 'viscosity='])

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
		elif option in ('-p', '--force'):
			force_profile = value
		elif option in ('-A', '--area'):
			A = float(value);	# in A^2 unit
		elif option in ('-v', '--viscosity'):
			v = float(value);	# in kg/m.s unit
	        elif option == NULL:
			print(usage)
			sys.exit(0)
	
	# more options: resname CNT, resname solvent, save BGF

	# main call
	getSlipLength(bgf_file, trj_file, force_profile)
