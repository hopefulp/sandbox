#!/home/noische/python

import sys, re, string, getopt, optparse, math, time, copy, pprint, os
from os import popen

import numpy as np

import bgf
import bgftools
import nutils as nu
import lammpstools as lt

usage = """
Write code description here. This is a help message.

filename.py -b bgf_file -t trj_file -o out_file -n #step 
"""
version = "150101"


def get_line(trj_file):
	"""
	Open a lazy file stream. This is considered to increase a file reading speed.
	"""
	with open(trj_file) as file:
		for i in file:
			yield i


def myfunction(bgf_file, trj_file, n_step, out_file, silent=False):

	# const
	PI = math.pi; vdw_r_C = 1.7;

	# init
	timestep = 0; l_timestep = []; line = []; n_header = 0;
	t1 = 0; t2 = 0; # clock

	# environments
	myBGF = bgf.BgfFile(bgf_file)	# load BGF
	myTRJ = open(trj_file)	# load LAMMPS trj
	myTRJ.seek(0)

	out_file = "myfunction.dat"
	f_out_file = open(out_file, 'w')	# data will be recorded here
	f_out_file.write(str(sys.argv) + "\n")

	curr_dir = os.path.abspath(".")
	temp_dir = curr_dir + "/data/"
	if not os.path.isdir(temp_dir): os.makedirs(temp_dir)
	print("* Current directory: " + curr_dir)
	print("* Temp directory: " + temp_dir)

	# variables
	pass;

	# scan LAMMPS trj to count how many shots in the file
	n_timestep = len(lt.getTrjInfo(trj_file))
	if n_step == 0:
		n_step = n_timestep;

	print("  ..The trajectory contains " + str(n_timestep) + " timesteps.")
	print("The script will proceed for the last " + str(n_step) + " timesteps.")

	# prepare something with BGF file
	for atom in myBGF.a:
		# do something
		pass;


	### Read trajectory file header
	while 1:
		templine = myTRJ.readline()
		line.append(templine.strip('\n').strip('ITEM: '))
		n_header += 1
		if "ITEM: ATOMS" in templine:
			break;

	# initial trajectory information
	timestep = int(line[1])
	natoms = int(line[3])
	boxsize = [line[5].split(' ')[0], line[5].split(' ')[1], line[6].split(' ')[0], line[6].split(' ')[1], line[7].split(' ')[0], line[7].split(' ')[1]]
	boxsize = [float(i) for i in boxsize]
	keywords = line[8].strip('ATOMS ')

	# launch LAMMPS trj file
	dumpatom = get_line(trj_file)


	### loop over trj and do the analysis
	processed_step = 0;	# counter
	t1 = t2 = 0; elapsed_time = 0;	# timer: calculate ETA

	while 1:
		# Show progress
		t1 = time.time();
		remaining_time = elapsed_time * (n_step - processed_step)
		sys.stdout.write('\r' + "Reading timestep.. " + str(timestep) + " (Remaining time: " + "{0:4.1f}".format(remaining_time) + " seconds, " + "{0:4.1f} minutes".format(remaining_time/60) + ")")
		sys.stdout.flush()

		# stop analysis if requested n_step reached
		if processed_step == n_step:
			break;

		# Read a shot from the trj file
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
				
		# pbc wrap if required
		myBGF = bgftools.periodicMoleculeSort(myBGF, 0, myBGF.CRYSTX)


		### do some analysis
		"""
		asdflkjhasdlkfjadslkfjasdlkvjbaldskjfha;lsiefhja;lsejkfblnkdsjvb
		"""


		### record data
		output = str(timestep) + '\t' + str(123) + '\n'	
		f_out_file.write(output)


		### update timer
		sys.stdout.flush()
		t2 = time.time()	# time mark
		elapsed_time = t2 - t1;
		processed_step += 1;

	print('')
	f_out_file.close()
	print("numbers of water molecules are written in " + out_file + " ..Done.")

	return 1

	### end of function


if __name__ == "__main__":
	bgf_file = ""; trj_file = ""; ff_file = ""; n_step = 0; 

	options, args = getopt.getopt(sys.argv[1:], 'hb:t:n:o:', ['help', 'bgf=', 'trj=', 'step=', 'out='])

	if len(sys.argv) < 2:
		print(version)
		print(usage)
		sys.exit(0)

	print "Requested options: " + str(options)

	for option, value in options:
	        if option in ('-h', '--help'):
	                print(usage)
			sys.exit(0)
	        elif option in ('-b', '--bgf'):
	                bgf_file = str(value)
	        elif option in ('-t', '--trj'):
	                trj_file = str(value)
		elif option in ('-n', '--step'):
			n_step = int(value)
		elif option in ('-o', '--out'):
			out_file = str(value)
	        elif option == NULL:
			print(usage)
			sys.exit(0)
	
	# main call
	myfunction(bgf_file, trj_file, out_file, n_step)
