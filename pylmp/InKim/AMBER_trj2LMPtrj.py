#!/home/noische/program/python27/bin/python

import sys
import re
import string
import getopt
import optparse
import time
import os
import re

# custom module
import nutils as nu

option = ""; args = ""; 
usage = """
getAMBERTrajectory.py: 

Usage: getAMBERTrajectory.py -c crd_file -v vel_file 

"""

version = "120730"

def sortkey(list):
	return list[0];

def get_line(trj_file):
	"""
	Open a lazy file stream. This is considered to increase a file reading speed.
	"""
	with open(trj_file) as file:
		for i in file:
			yield i

def getAMBERTrajectory(crd_file, vel_file, LMPtrj_file, silent=False):
	"""
	Briefings:
		1. read mdcrd in a big list: omit volume info
		2. read mdvel in a big list
		3. write into lammpstrj atom by atom
	"""
	### initialize
	boxsize = []; cntLine = 0; coord = []; natoms = 0;
	temp_crd_file = ".temp.mdcrd"
	temp_vel_file = ".temp.mdvel"

	### Preprocess mdcrd: delete the title and save it into a temp file
	if os.path.getsize(crd_file) == 0: nu.die(crd_file + " is empty.")
	if not silent: print("Preprocessing " + temp_crd_file + " .. this may take a while.")
	temp_crd_cmd = "sed '1d' " + crd_file + " > " + temp_crd_file
	os.system(temp_crd_cmd)

	### Preprocess mdvel: delete the title and save it into a temp file
	if os.path.getsize(vel_file) == 0: nu.die(vel_file + " is empty.")
	if not silent: print("Preprocessing " + temp_vel_file + " .. this may take a while.")
	temp_vel_cmd = "sed '1d' " + vel_file + " > " + temp_vel_file
	os.system(temp_vel_cmd)

	### find file length of one step in mdcrd
	f_crd_file = open(temp_crd_file)
	try:
		line = f_crd_file.readline()	# blow up the 1st line: title
	except:
		nu.die("Check the mdcrd file " + str(crd_file))

	while 1:
		line = f_crd_file.readline()
		if not line: break	

		pattern = re.compile("\s*(\d*.\d*)\s*(\d*.\d*)\s*(\d*.\d*)$")
		cntLine += 1;
		if re.match(pattern, line):
			cntLine += 1;
			break;
		
	f_crd_file.close()

	# test
	#print("test: cntLine " + str(cntLine))

	### load coordinate and velocity
	# REMARK: every item in coord exactly matches to a shot, whereas every item in vel doesn't.
	if not silent: print("Loading coordinates and velocity file.")
	coord = nu.openBigFile(temp_crd_file, cntLine, True, False)	
	vel = nu.openBigFile(temp_vel_file, cntLine - 1, True, False)	# might be a one shot.

	if not silent: print("Coordinates and velocities are loaded.")

	### get box information from coord
	for i in coord:
		boxsize.append(i.pop())

	### check the number of boxsize and shots in coord
	if len(boxsize) != len(coord):
		nu.die("Failed to read boxsizes in the mdcrd file.")

	#test# print("len(coord[0]) "+str(len(coord[0])))	# check

	### how many atoms in there?
	nitems_crd = 0; nitems_vel = 0; natoms_crd = 0; natoms_vel = 0; nitems = 0; 
	temp_crd = [ item for sublist in coord[0] for item in sublist ]
	nitems_crd = len(temp_crd)
	if nitems_crd % 3 == 0:
		natoms_crd = nitems_crd / 3
	else:
		nu.die("The number of atoms in " + crd_file + " is suspicious. Please check the file.")
	del(temp_crd)

	# test
	print("test:natoms_crd: " + str(natoms_crd))

	temp_vel = [ item for sublist in vel[0] for item in sublist ]
	nitems_vel = len(temp_vel)
	if nitems_vel % 3 == 0:
		natoms_vel = nitems_vel / 3
	else:
		nu.die("The number of atoms in " + vel_file + " is suspicious. Please check the file.")
	del(temp_vel)

	# test
	print("test:natoms_vel: " + str(natoms_vel))

	if natoms_crd == natoms_vel:
		natoms = natoms_crd
	else:
		nu.die("Atom number mismatch in " + crd_file + " and " + vel_file + ".")

	# test
	#print(nitems)

	### make a lammpstrj format file
	f_LMPtrj_file = open(LMPtrj_file, 'w')
	for index, shot in enumerate(coord):
		output = "";

		# timestep = index of the shot in coord
		output += "ITEM: TIMESTEP" + "\n"
		output += str(index) + "\n"

		# number of atoms = natoms_crd
		output += "ITEM: NUMBER OF ATOMS" + "\n"
		output += str(natoms) + "\n"

		# box bounds = boxsize
		output += "ITEM: BOX BOUNDS" + "\n"
		output += "0 " + str(boxsize[index][0]) + "\n"
		output += "0 " + str(boxsize[index][1]) + "\n"
		output += "0 " + str(boxsize[index][2]) + "\n"

		# atoms
		output += "ITEM: ATOMS id type xu yu zu vx vy vz" + "\n"


	## end of getLAMMPSTrajectory


if __name__ == "__main__":

	print("\n" + sys.argv[0] + " version " + str(version) + "\n")

	if len(sys.argv) < 2:
		print(usage);
		sys.exit(1)

	crd_file = ""; vel_file = ""; LMPtrj_file = "";
	timestep = -1;

	options, args = getopt.getopt(sys.argv[1:], 'hc:v:t:', ['help','crd=','vel=','trj='])
	try:
		for option, value in options:
			if option in ('-h', '--help'):
				print(usage)
				sys.exit(0);
			elif option in ('-c', '--crd'):
				crd_file = value
			elif option in ('-v', '--vel'):
				vel_file = value
			elif option in ('-t', '--trj'):
				LMPtrj_file = value
			elif option == NULL:
				print(usage)
				sys.exit(0)
	except ValueError:
		nu.die("Script cannot continue. Check input parameters.")

	getAMBERTrajectory(crd_file, vel_file, LMPtrj_file, silent=False)


