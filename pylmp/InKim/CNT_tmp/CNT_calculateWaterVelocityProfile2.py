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

option = ""; args = ""; bgf_file = ""; trj_file = ""; out_file = ""; n_skip = 0;
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

def getSlipLength(bgf_file, trj_file, ff_file, interval, n_step, n_avg_step, silent=False):
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

	myBGF = bgf.BgfFile(bgf_file)
	myTRJ = open(trj_file)
	myTRJ.seek(0)
	myPkl = open(trj_file + ".profile.r_ave_vel-s" + str(n_step) + "k" + str(n_skip) + ".pickle", 'w')
	#f_dump = open(trj_file + ".profile.r_ave_vel.dump", 'w')
	#dump = "";

	# how many steps to go?
	n_timestep = len(lt.getTrjInfo(trj_file, silent=False))
	if n_step == 0:
		n_step = n_timestep;

	if not silent: print("The trajectory contains " + str(n_timestep) + " timesteps.")

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

	# calculate radius_CNT
	avg_radius_CNT = 0.0; x_CNT = 0.0; y_CNT = 0.0;
	if not silent: print("Calculating the average radius of CNT throughout the averaging steps")
	while 1:

		### Show progress
		t1 = time.time();
		remaining_time = elapsed_time * (n_step - processed_step)

		### Read
		try:
			chunk = [next(dumpatom) for i in range(natoms+n_header)]
		except StopIteration:
			break;

		# skip if necessary
		if processed_step <= n_skip:
			processed_step += 1;
			continue;

		timestep = int(chunk[1])
		l_timestep.append(timestep)
		natoms = int(chunk[3])
		boxsize = [chunk[5].split(' ')[0], chunk[5].split(' ')[1], chunk[6].split(' ')[0], chunk[6].split(' ')[1], chunk[7].split(' ')[0], chunk[7].split(' ')[1]]; 
		boxsize = [float(i) for i in boxsize]; boxsize = [(boxsize[1] - boxsize[0]), (boxsize[3] - boxsize[2]), (boxsize[5] - boxsize[4])]
		keywords = chunk[8].split('ATOMS ')[1].strip('\n').split(' ')
		sys.stdout.write('\r' + "Fetched timestep: " + str(timestep) + " (Remaining time: " + "{0:4.1f}".format(remaining_time) + " seconds, " + "{0:4.1f} minutes".format(remaining_time/60) + ")")
		sys.stdout.flush()

		### update myBGF with trajectory information ###
		natom_bgf = len(myBGF.a)	# number of atoms in BGF file
	
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

			try:
				for i in range(0, 3):
					myBGF.CRYSTX[i] = boxsize[i]
			except:
				pass;
		
		# apply periodic condition
		#myBGF = bgftools.periodicMoleculeSort(myBGF, 0, True)

		### myBGF update complete! ###

		# for CNT atoms, calculate some properties
		min_x_CNT = 0.0; max_x_CNT = 0.0; radius_CNT = 0.0; height_CNT = 0.0; aNo_WAT_O_not_in_CNT = []; temp = [];
		min_y_CNT = 0.0; max_y_CNT = 0.0;
		min_z_CNT = 0.0; max_z_CNT = 0.0;
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

		if processed_step > 3:
			break;

	avg_radius_CNT, s = nu.meanstdv(l_avg_radius_CNT)
	bins = np.arange(0.0, math.ceil(radius_CNT), interval);
	l_radial_vz_profile = np.zeros(len(bins)-1)

	if not silent: print("\nAveraged CNT radius: " + str(avg_radius_CNT) + " A")
	if not silent: print("CNT center: " + str([x_CNT, y_CNT]))


	# calculate properties
	myTRJ.close()
	myTRJ = open(trj_file)
	myTRJ.seek(0)
	dumpatom = get_line(trj_file)
	processed_step = 0;
	r = []; vz = [];	# data container
	while 1:
		data = dict();	# stores [r vz] for a timestep

		### Show progress
		t1 = time.time();
		remaining_time = elapsed_time * (n_step - processed_step)

		### Read
		try:
			chunk = [next(dumpatom) for i in range(natoms+n_header)]
		except StopIteration:
			break;

		if processed_step <= n_skip:
			processed_step += 1;
			continue;

		timestep = int(chunk[1])
		l_timestep.append(timestep)
		natoms = int(chunk[3])
		boxsize = [chunk[5].split(' ')[0], chunk[5].split(' ')[1], chunk[6].split(' ')[0], chunk[6].split(' ')[1], chunk[7].split(' ')[0], chunk[7].split(' ')[1]]; 
		boxsize = [float(i) for i in boxsize]; boxsize = [(boxsize[1] - boxsize[0]), (boxsize[3] - boxsize[2]), (boxsize[5] - boxsize[4])]
		keywords = chunk[8].split('ATOMS ')[1].strip('\n').split(' ')
		sys.stdout.write('\r' + "Fetched timestep: " + str(timestep) + " (Remaining time: " + "{0:4.1f}".format(remaining_time) + " seconds, " + "{0:4.1f} minutes".format(remaining_time/60) + ")")
		sys.stdout.flush()

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


		# calculate z-directional velocity of water molecules wrt r		vcm = sum vi mi / sum mi
		atom = myBGF.getAtom(aNo_WAT_O[0])
		mass_O = bgftools.getMass(myBGF, [aNo_WAT_O[0]], ff_file)
		atom = myBGF.getAtom(atom.CONECT[0])
		mass_H = bgftools.getMass(myBGF, [atom.aNo], ff_file)
		mass_H2O = mass_O + 2 * mass_H	# H2O mass


		for ano in aNo_WAT_O:
			atom = myBGF.getAtom(ano)	# oxygen

			# com & velocity calculation
			cmx = atom.x * mass_O
			cmy = atom.y * mass_O
			vcmz = atom.vz * mass_O

			for ano2 in atom.CONECT:
				atom2 = myBGF.getAtom(ano2)
				cmx += atom2.x * mass_H
				cmy += atom2.y * mass_H
				vcmz += atom2.vz * mass_H

			cmx /= mass_H2O
			cmy /= mass_H2O
			vcmz /= mass_H2O

			dist = math.sqrt((x_CNT - cmx)**2 + (y_CNT - cmy)**2)
			r.append(dist)
			vz.append(vcmz*10**5)	# LAMMPS real unit conversion for velocity
			

		# timestep mark for starting average
		if l_avg_count == []:
			n_avg_count_step = timestep

		l_avg_count.append(timestep)


		# if the time has come.. average!
		if len(l_avg_count) == n_avg_step or processed_step == n_timestep - 1:
			if not silent: print(" Averaging invoked at timestep " + str(timestep) + " (" + str(len(l_avg_count)) + " points)")

			# binning
			vz = np.array(vz)
			sum_vz = np.ma.array(np.histogram(r, bins, weights=vz)[0])
			sum_vz2 = np.ma.array(np.histogram(r, bins, weights=vz*vz)[0])
			n_vz = np.ma.array(np.histogram(r, bins)[0])
			#mask_n = ma.masked_values(n_vz, 0)

			mean = sum_vz / n_vz
			std = np.sqrt(sum_vz2/n_vz - mean*mean)

			mean = np.ma.fix_invalid(mean, fill_value=0.0).data
			std = np.ma.fix_invalid(std, fill_value=0.0).data
			n_vz = np.ma.fix_invalid(n_vz, fill_value=0).data

			if len(mean) != len(bins)-1 or len(std) != len(bins)-1:
				nu.die("on numpy masked_array calculation!")

			l_result.append([l_avg_count, bins, n_vz, mean, std])

			# reset
			r = []; vz = [];
			l_radial_vz_profile = np.zeros(len(bins)-1)
			l_avg_count = [];


		t2 = time.time()	# time mark
		elapsed_time = t2 - t1;
		processed_step += 1

	#f_dump.write(dump)
	pkl.dump(l_result, myPkl)

	print('')
	return 1

	### end of function


if __name__ == "__main__":
	# n_avg_step = how many steps will be averaged
	bgf_file = ""; trj_file = ""; ff_file = ""; n_step = 0; n_avg_step = 1; interval = 0.2; n_pass = 0;

	options, args = getopt.getopt(sys.argv[1:], 'hb:t:n:s:f:i:k:', ['help', 'bgf=', 'trj=', 'number=', 'step=', 'ff=', 'interval=', 'skip='])

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
		elif option in ('-i', '--interval'):
			interval = float(value)
		elif option in ('-k', '--skip'):
			n_skip = int(value)
	        elif option == NULL:
			print(usage)
			sys.exit(0)
	
	# more options: resname CNT, resname solvent, save BGF

	# main call
	getSlipLength(bgf_file, trj_file, ff_file, interval, n_step, n_avg_step)
