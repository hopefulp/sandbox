#!/home/noische/program/python27/bin/python

import sys
import re
import string
import getopt
import optparse
import time
import os
import math

# BGF module
import bgf
import bgftools
import nutils as nu

# numpy
import numpy

option = ""; args = ""; 
usage = """
checkMETalign.py: 

Usage: checkMETalign.py -b bgf_file -t trj_file -o out_file -n timestep -p (periodic sort)

"""

version = "120911"
#trj_file = ""
def sortkey(list):
	return list[0];

def checkMETalign(bgf_file, out_file, verbose, silent=False):
	"""
	This is the script only for checking the MET alignment status in (6,6) CNT.
	Written from 2012. 09. 11.
	Only considers for methanol_CNT cases.
	"""

	timestep = 0;
	natoms = 0;
	boxsize = [0, 0, 0, 0, 0, 0];
	keywords = [];
	data = [];
	line = []
	l_timestep_seek = []
	n_header = 0;
	mode = "";

	# Open BGF
	if isinstance(bgf_file, bgf.BgfFile):
		myBGF = bgf_file
	else:
		if not silent: print("Reading " + bgf_file + " ..")
		myBGF = bgf.BgfFile(bgf_file)

	# extract aNos of CNT and MET
	aNo_CNT = []; aNo_MtOH_C = []
	for atom in myBGF.a:
		if "CNT" in atom.rName:
			aNo_CNT.append(atom.aNo)
		if "MET" in atom.rName and "C" in atom.aName:
			aNo_MtOH_C.append(atom.aNo)

	N_CNT = len(aNo_CNT)

	# find orientation of CNT
	min_x_CNT = 0.0; max_x_CNT = 0.0; 
        min_y_CNT = 0.0; max_y_CNT = 0.0;
        min_z_CNT = 0.0; max_z_CNT = 0.0;
	radius_CNT = 0.0; height_CNT = 0.0; aNo_MtOH_C_not_in_CNT = []; temp = [];
        CNT_orientation = "";

	for aNo in aNo_CNT:
		atom = myBGF.getAtom(aNo)

		if atom.x < min_x_CNT:
			min_x_CNT = atom.x
		if atom.x > max_x_CNT:
			max_x_CNT = atom.x
		if atom.y < min_y_CNT:
			min_y_CNT = atom.y
		if atom.y > max_y_CNT:
			max_y_CNT = atom.y
		if atom.z < min_z_CNT:
			min_z_CNT = atom.z
		if atom.z > max_z_CNT:
			max_z_CNT = atom.z
	
	x_diff = max_x_CNT - min_x_CNT
	y_diff = max_y_CNT - min_y_CNT
	z_diff = max_z_CNT - min_z_CNT

	if x_diff > y_diff and x_diff > z_diff:
		# CNT is aligned along x axis
		height_CNT = x_diff
		CNT_orientation = "x"
	elif y_diff > x_diff and y_diff > z_diff:
		# CNT is aligned along y axis
		height_CNT = y_diff
		CNT_orientation = "y"
	elif z_diff > x_diff and z_diff > y_diff:
		# CNT is aligned along z axis
		height_CNT = z_diff
		CNT_orientation = "z"

	if CNT_orientation == "x":
		for aNo in aNo_CNT:
			atom = myBGF.getAtom(aNo)
			radius_CNT += math.sqrt(atom.y**2 + atom.z**2)	# average radius of CNT
	elif CNT_orientation == "y":
		for aNo in aNo_CNT:
			atom = myBGF.getAtom(aNo)
			radius_CNT += math.sqrt(atom.x**2 + atom.z**2)	# average radius of CNT
	elif CNT_orientation == "z":
		for aNo in aNo_CNT:
			atom = myBGF.getAtom(aNo)
			radius_CNT += math.sqrt(atom.x**2 + atom.y**2)	# average radius of CNT

	radius_CNT = radius_CNT / N_CNT

	### Rotate y or z to x axis
	Tyx = numpy.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
	Tzx = numpy.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])

	if CNT_orientation == "y":
		for atom in myBGF.a:
			u = numpy.matrix([atom.x, atom.y, atom.z]).T
			Tu = Tyx*u
			atom.x = float(Tu[0])
			atom.y = float(Tu[1])
			atom.z = float(Tu[2])

	elif CNT_orientation == "z":
		for atom in myBGF.a:
			u = numpy.matrix([atom.x, atom.y, atom.z]).T
			Tu = Tzx*u
			atom.x = float(Tu[0])
			atom.y = float(Tu[1])
			atom.z = float(Tu[2])

	### check MET position

	output = "";
	MeOH_align = [];

	# find O in MET: MeOH_align = [C.x, O.x]
	for ano in aNo_MtOH_C:
		atom = myBGF.getAtom(ano)
		for i in atom.CONECT:
			atom2 = myBGF.getAtom(i)
			if "O" in atom2.ffType:
				MeOH_align.append([atom.x, atom2.x])

	MeOH_align.sort(key=sortkey)

	for i in MeOH_align:
		if i[0] > i[1]:
			# C-O: 1
			if verbose:
				output += "{0:-8.5f}".format(i[0]) + "\t" + "{0:-8.5f}".format(i[1]) + "\t" + "1\n"
			else:
				output += "1 "
		else:
			# O-C: 0
			if verbose:
				output += "{0:-8.5f}".format(i[0]) + "\t" + "{0:-8.5f}".format(i[1]) + "\t" + "0\n"
			else:
				output += "0 "

	output += "\n"

	print(output)
		
	# save
	#myBGF.saveBGF(out_file)


	## end of checkMETalign


if __name__ == "__main__":

	print("\n" + sys.argv[0] + " version " + str(version) + "\n")

	if len(sys.argv) < 2:
		print(usage);
		sys.exit(1)

	bgf_file = ""; trj_file = ""; out_file = ""; verbose = False;
	ano = -1;

	options, args = getopt.getopt(sys.argv[1:], 'hb:o:v:', ['help','bgf=','out=', 'verbose='])
	try:
		for option, value in options:
			if option in ('-h', '--help'):
				print(usage)
				sys.exit(0);
			elif option in ('-b', '--bgf'):
				bgf_file = value
			elif option in ('-o', '--out'):
				out_file = value
			elif option in ('-v', '--verbose'):
				verbose = true
			elif option == NULL:
				print(usage)
				sys.exit(0)
	except ValueError:
		nu.die("script cannot continue. Check input parameters.")

	if out_file == "":
		out_file = bgf_file.split(".bgf")[0] + ".check.bgf"

	checkMETalign(bgf_file, out_file, verbose, silent=False)


