#!/home/noische/python

import sys, re, string, getopt, optparse, math, time, pprint, os
from os import popen
import bgf
import bgftools
import numpy
import operator
import copy
import nutils as nu
import lammpstools as lt


option = ""; args = ""; dat_file = ""; trj_file = ""; out_file = ""; atom_type = ""; criteria_distance = 0;
usage = """
LAMMPS_trj2data.py -b bgf_file -d dat_file -t trj_file -s step (optional to get bgf)
"""
version = "130614"

#_________________
# 130603: first version
# 130614: applies box information from trajectory
#_________________

def get_line(trj_file):
	"""
	Open a lazy file stream. This is considered to increase a file reading speed.
	"""
	with open(trj_file) as file:
		for i in file:
			yield i

def trj2data(dat_file, trj_file, step, skipVelocity, silent=False):

	# init
	timestep = 0; l_timestep = []; line = []; n_header = 0;
	t1 = 0; t2 = 0; # clock

	myBGF = bgf.BgfFile(bgf_file)
	myTRJ = open(trj_file)
	myTRJ.seek(0)

	### LAMMPS Data
	dat_atoms = [];
	if not silent: print("== Step 1. Loading LAMMPS Data")

	if not os.path.exists(dat_file):
		nu.die("Please check the LAMMPS data file.")

	f_dat_file = open(dat_file)
	temp = f_dat_file.read().split('\n')
	dat_atom_start = temp.index('Atoms')
	dat_atom_end = temp.index('Bonds')
	dat_atoms = temp[dat_atom_start + 1:dat_atom_end]
	f_dat_file.close()


	### LAMMPS Trajectory
	# how many steps to go?
	if not silent: print("== Step 2. Loading LAMMPS Trajectory")
	n_timestep = len(lt.getTrjInfo(trj_file))
	print("The trajectory contains " + str(n_timestep) + " timesteps.")

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

	if not silent: print("== Step 3. Reading LAMMPS Trajectory")
	while 1:
		try:
			chunk = [next(dumpatom) for i in range(natoms+n_header)]
		except StopIteration:
			break;
	
		timestep = int(chunk[1])
		natoms = int(chunk[3])
		boxsize = [chunk[5].split(' ')[0], chunk[5].split(' ')[1], chunk[6].split(' ')[0], chunk[6].split(' ')[1], chunk[7].split(' ')[0], chunk[7].split(' ')[1]]; 
		boxsize = [float(i) for i in boxsize]; boxsize = [(boxsize[1] - boxsize[0]), (boxsize[3] - boxsize[2]), (boxsize[5] - boxsize[4])]
		keywords = chunk[8].split('ATOMS ')[1].strip('\n').split(' ')

		### Show progress
		t1 = time.time();
		remaining_time = elapsed_time * (n_timestep - processed_step)
		sys.stdout.write('\r' + "Reading timestep.. " + str(timestep) + " (Elapsed time for the previous step: " + "{0:4.1f}".format(elapsed_time) + " seconds, Remaining time: " + "{0:4.1f}".format(remaining_time) + " seconds = " + "{0:4.1f} minutes".format(remaining_time/60) + ")")
		sys.stdout.flush()

		processed_step += 1;

		### if step is specified, then save only for that step.
		if step != "":
			if timestep == step:
				pass;
			else:
				continue;
		else:
			pass;

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

		flagVelocity = False;
		if 'vx' in keywords or 'vy' in keywords or 'vz' in keywords:
			flagVelocity = True
	
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

			if not skipVelocity:
				if flagVelocity == True:
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
		myBGF = bgftools.periodicMoleculeSort(myBGF, 0)

		#myBGF.saveBGF(bgf_file[:-4] + "." + str(timestep) + ".trjupdated.bgf")
		### myBGF update complete! ###

		### update LAMMPS data file ###
		myDAT = open(dat_file)
		myNewDAT = open(dat_file + "." + str(timestep), 'w')
		linecount = 0;

		while 1:
			line = myDAT.readline()

			if not line:
				break;

			if linecount == 0:
				output = line.strip("\n") + " and updated with timestep " + str(timestep) + " in " + trj_file + "\n"
				linecount += 1;
				myNewDAT.write(output)
				continue;

			if "xlo " in line:
				output = "\t" + " {0:>10.6f} {1:>10.6f}".format(0.0, boxsize[0]) + " xlo xhi\n"
				myNewDAT.write(output)
				continue;
			elif "ylo " in line:
				output = "\t" + " {0:>10.6f} {1:>10.6f}".format(0.0, boxsize[1]) + " ylo yhi\n"
				myNewDAT.write(output)
				continue;
			elif "zlo " in line:
				output = "\t" + " {0:>10.6f} {1:>10.6f}".format(0.0, boxsize[2]) + " zlo zhi\n"
				myNewDAT.write(output)
				continue;

			if not "Atoms" in line:
				myNewDAT.write(line)

			else:
				myNewDAT.write("Atoms\n\n")

				for i in dat_atoms:
					if len(i) == 0:
						continue;

					parse = i.split()
					atomNo = int(parse[0])
					molNo = int(parse[1])
					atomTypeNo = int(parse[2])

					atom = myBGF.getAtom(atomNo)
					if atom.aNo != atomNo or atom.rNo != molNo:
						nu.die("BGF and data file mismatch. Please check both files.")
					else:
						output = "{0:>8} {1:>8} {2:>8}  {3:10.5f} {4:10.5f} {5:10.5f} {6:10.5f} {7:12.8f} {8:12.8f} {9:12.8f}".format(atom.aNo, atom.rNo, atomTypeNo, atom.charge, atom.x, atom.y, atom.z, atom.vx, atom.vy, atom.vz) + "\n"
					myNewDAT.write(output)

				# consume
				for i in range(len(dat_atoms)):
					line = myDAT.readline()

				myNewDAT.write("\n")


		t2 = time.time()	# time mark
		elapsed_time = t2 - t1;

	print('')
	return 1

	### end of function


if __name__ == "__main__":
	bgf_file = ""; dat_file = ""; trj_file = ""; step = ""; skipVelocity = False;

	options, args = getopt.getopt(sys.argv[1:], 'hb:d:t:s:v', ['help', 'bgf=', 'dat=', 'trj=', 'step=', 'vel'])

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
	        elif option in ('-d', '--dat'):
	                dat_file = value
	        elif option in ('-t', '--trj'):
	                trj_file = value
		elif option in ('-s', '--step'):
			step = int(value)
		elif option in ('-v', '--vel'):
			skipVelocity = True
	        elif option == NULL:
			print(usage)
			sys.exit(0)
	
	# more options: resname CNT, resname solvent, save BGF

	# main call
	trj2data(dat_file, trj_file, step, skipVelocity)
