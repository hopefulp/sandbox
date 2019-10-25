#!/opt/applic/epd/bin/python

import sys
import re
import string
import getopt
import optparse
import math
import time

from os import popen
import bgf

#-----------------
# update the coordinate in the original BGF file from LAMMPS trajectory file
#_________________
def createPeriodicBox(bgf_file, out_file, silent=True):
	boxsize = [0, 0, 0, 0, 0, 0]

	# open bgf
	if isinstance(bgf_file, bgf.BgfFile):
		myBGF = bgf_file
	else:
		if not silent: print("reading " + bgf_file + " ..")
		myBGF = bgf.BgfFile(bgf_file)

	xlo = ylo = zlo = 1000; xhi = yhi = zhi = -1000;

	for atom in myBGF.a:
		if "C_2G" in atom.ffType or "CA" in atom.ffType:
			if atom.x < xlo:
				xlo = atom.x
			if atom.x > xhi:
				xhi = atom.x
			if atom.y < ylo:
				ylo = atom.y
			if atom.y > yhi:
				yhi = atom.y

		if atom.z < zlo:
			zlo = atom.z
		if atom.z > zhi:
			zhi = atom.z
		
	"""
	for atom in myBGF.a:
		if atom.x < boxsize[0]:
			boxsize[0] = atom.x
		elif atom.x > boxsize[1]:
			boxsize[1] = atom.x

		if atom.y < boxsize[2]:
			boxsize[2] = atom.y
		elif atom.y > boxsize[3]:
			boxsize[3] = atom.y

		if atom.z < boxsize[4]:
			boxsize[4] = atom.z
		elif atom.z > boxsize[5]:
			boxsize[5] = atom.z
	"""

	### Size of the box
	a = []
	"""
	a.append(boxsize[1] - boxsize[0] + 1)	# 4 for vdw
	a.append(boxsize[3] - boxsize[2] + 1)	# "
	a.append(boxsize[5] - boxsize[4] + 1)	# "
	"""
	a.append(xhi-xlo+1)
	a.append(yhi-ylo+1)
	a.append(zhi-zlo+1)
	a.append(90.0)
	a.append(90.0)
	a.append(90.0)

	myBGF.PERIOD = "111"
	myBGF.AXES = "ZYX"
	myBGF.SGNAME = "P 1                  1    1\n"
	myBGF.CELLS = [-1, 1, -1, 1, -1, 1]
	myBGF.CRYSTX = a

	# save
	if isinstance(out_file, str):
		if not silent: print("saving information to " + out_file + " ..")
		myBGF.saveBGF(out_file)
		return 1;
	else:
		return myBGF;


if __name__ == "__main__":

	option = ""; args = ""; bgf_file = ""; out_file = "";
	usage = """
Usage: createPeriodicBox.py -b bgf_file -o out_file
	Creates a periodic box to a BGF file.
	It is highly recommended to center your BGF with ~tpascal/scripts/centerbgf.pl
	"""

        if len(sys.argv) < 2:
                print(usage)
		sys.exit(0)

	options, args = getopt.getopt(sys.argv[1:], 'hb:o:', ['help','bgf=','out='])
	for option, value in options:
	        if option in ('-h', '--help'):
	                print(usage); sys.exit(0)
	        elif option in ('-b', '--bgf'):
	                bgf_file = value
		elif option in ('-o', '--out'):
			out_file = value
	        elif option in NULL:
	                print(usage); sys.exit(0)

	# main call
	createPeriodicBox(bgf_file, out_file)
