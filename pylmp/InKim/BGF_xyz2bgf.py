#!/opt/applic/epd/bin/python

import sys
import re
import string
import getopt
import optparse
import math
import os

import bgf
import nutils as nu

version = "130604"

def xyz2bgf(xyz_file, bgf_file, silent=False):
	"""
	Converts xyz to bgf
	"""

	### open xyz
	f_xyz = open(xyz_file, 'r')

	### create bgf
	myBGF = bgf.BgfFile()
	myBGF.BIOGRF = '200'
	myBGF.DESCRP = xyz_file

	### read and create atoms for bgf file
	while 1:
		line = f_xyz.readline()

		if not line: break;

		words = string.split(line)

		if not words: break;

		natoms = int(words[0])
		myBGF.REMARK = [cleansym(f_xyz.readline())]

		for i in range(natoms):
			line = f_xyz.readline()
			words = string.split(line)
			sym = words[0]
			sym = cleansym(sym)
			x, y, z = float(words[1]), float(words[2]), float(words[3])

			### add atoms
			atom = bgf.BgfAtom()
			atom.x = x; 
			atom.y = y; 
			atom.z = z; 
			atom.ffType = sym;
			atom.aNo = i;
			atom.aName = sym + str(i)
			atom.rNo = 1;
			atom.chain = "A"
			atom.rName = "RES"
			myBGF.addAtom(atom)

	### save bgf
	myBGF.OTHER = []
	myBGF.saveBGF(bgf_file)

	### end of function


def cleansym(s):
	return re.split('[^a-zA-Z]',s)[0]


if __name__ == "__main__":
	bgf_file = ""; qchem_in_file = ""; rem_file = ""; 
	usage = """
BGF_xyz2bgf.py -i xyz_file -o bgf_file
	"""

	options, args = getopt.getopt(sys.argv[1:], 'hi:o:', ['help', 'xyz=', 'bgf='])

	if len(sys.argv) < 2:
		print(usage)
		sys.exit(0)

	print("Requested options: " + str(options))

	for option, value in options:
	        if option in ('-h', '--help'):
	                print(usage)
			sys.exit(0)
	        elif option in ('-o', '--bgf'):
	                bgf_file = value
	        elif option in ('-i', '--xyz'):
	                xyz_file = value
	        elif option == NULL:
			print(usage)
			sys.exit(0)

	# required options
	if xyz_file == "" or os.path.exists(xyz_file) == False:
		nu.die("No XYZ file named " + xyz_file)

	# default options
	if bgf_file == "":
		bgf_file = xyz_file.split(".xyz")[0] + ".in"
	
	# main call
	xyz2bgf(xyz_file, bgf_file)
