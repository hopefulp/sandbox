#!/home/noische/python

import sys, re, string, getopt
import os
import bgf
import math
import numpy

option = ""; args = ""; bgf_file = ""; trj_file = ""; out_file = ""; atom_type = ""; criteria_distance = 0;
usage = """
CNT_adjustSimBoxSize.py -b bgf_file -f ff_file
"""

version = "150204"

#-----------------
# adjust the simulation box size to counteract the PPPM problem during NPT
# which is considered due to the long bond length
# by in kim 
# last update: 2012/06/27
#_________________
def adjustSimboxSize(bgf_file, ff_file, rotate, silent=False):

	# init
	radius_CNT = 0.0; height_CNT = 0.0; temp = [];
	min_x_CNT = 10000.0; max_x_CNT = -10000.0; 
	min_y_CNT = 10000.0; max_y_CNT = -10000.0;
	min_z_CNT = 10000.0; max_z_CNT = -10000.0;
	CNT_orientation = "";

	### check whether CNT or BNNT
	myBGF = bgf.BgfFile(bgf_file)
	resNameCNT = ""
	for atom in myBGF.a:
		if "CNT" in atom.rName:
			resNameCNT = "CNT"
			break;
		elif "BNT" in atom.rName:
			resNameCNT = "BNT"

	### basic ff file
	if resNameCNT == "CNT":
		ff_file = "/home/noische/ff/graphite_lj.ff " + ff_file
	else:
		ff_file = "/home/noische/ff/DREIDING2.21.ff " + ff_file


	# extract aNos of CNT in the BGF file
	aNo_CNT = []
	for atom in myBGF.a:
		if resNameCNT in atom.rName:
			aNo_CNT.append(atom.aNo)

	N_CNT = len(aNo_CNT)	# the number of CNT atoms

	# for CNT atoms, calculate some properties
	# check the orientation of CNT
	for aNo in aNo_CNT:
		atom = myBGF.getAtom(aNo)
		# minimum and maximum x coord of CNT: this will be the height of CNT
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

	if rotate:
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
				min_x_CNT = min_y_CNT;
				max_x_CNT = max_y_CNT;

		elif CNT_orientation == "z":
			for atom in myBGF.a:
				u = numpy.matrix([atom.x, atom.y, atom.z]).T
				Tu = Tzx*u
				atom.x = float(Tu[0])
				atom.y = float(Tu[1])
				atom.z = float(Tu[2])
				min_x_CNT = min_z_CNT;
				max_x_CNT = max_z_CNT;

		# CNT is aligned along x axis
		myBGF.CRYSTX = [height_CNT + 0.2, radius_CNT + 30.0, radius_CNT + 30.0, 90.0, 90.0, 90.0];
	else:
		if CNT_orientation == "x":
			myBGF.CRYSTX = [height_CNT + 0.2, radius_CNT + 30.0, radius_CNT + 30.0, 90.0, 90.0, 90.0];
		elif CNT_orientation == "y":
			myBGF.CRYSTX = [radius_CNT + 30.0, height_CNT + 0.2, radius_CNT + 30.0, 90.0, 90.0, 90.0];
		elif CNT_orientation == "z":
			myBGF.CRYSTX = [radius_CNT + 30.0, radius_CNT + 30.0, height_CNT + 0.2, 90.0, 90.0, 90.0];

	myBGF.PERIOD = "111"
	myBGF.AXES = "ZYX"
	myBGF.SGNAME = "P 1                  1    1"
	myBGF.CELLS = [-1, 1, -1, 1, -1, 1]

	# write output
	save_filename = bgf_file[:-4] + ".margined.bgf"
	myBGF.saveBGF(save_filename)
	centerBGF_cmd = "~tpascal/scripts/centerBGF.pl -b " + save_filename + " -f '" + ff_file + "' -c com_center "
	os.system(centerBGF_cmd)

	return 1

	### end of function


if __name__ == "__main__":
	bgf_file = ""; ff_file = ""; rotate = False;

	options, args = getopt.getopt(sys.argv[1:], 'hb:f:r', ['help', 'bgf=', 'ff=', 'rotate='])

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
		elif option in ('-f', '--ff'):
			ff_file = value
		elif option in ('-r', '--rotate'):
			rotate = True
	        elif option == NULL:
			print(usage)
			sys.exit(0)
	
	# main call
	adjustSimboxSize(bgf_file, ff_file, rotate)
