#!/opt/applic/epd/bin/python
"""
PEI_NNstrain.py
Original: Aug 23 2011 In Kim
"""

# Python Modules
import sys
import os
import string
import getopt
import time
import math

# Custom Modules
import numpy as np
import bgf
import bgftools

# Globals
version = '110823'


def rotate(bgf_file, v, out_file, silent=True):
	"""
	rotate BGF to v2 to z axis.
	"""
	# Initialization

	# Open BGF
	if isinstance(bgf_file, bgf.BgfFile):
		myBGF = bgf_file
	else:
		if not silent: print("opening bgf file.. " + str(bgf_file))
		myBGF = bgf.BgfFile(bgf_file)

	# rotation 1 to x axis
	v1_dot_v2 = v[0]
	mag_v1 = 1.0
	mag_v2 = math.sqrt(v[0]**2 + v[1]**2)
	theta = math.acos(v1_dot_v2 / (mag_v1 * mag_v2))
	theta = -theta
	print(theta)
	rot1 = np.array([[math.cos(theta), -math.sin(theta), 0], [math.sin(theta), math.cos(theta), 0], [0, 0, 1]])

	# rotation 2 to y axis
	v1_dot_v2 = v[1]
	mag_v1 = 1.0
	mag_v2 = math.sqrt(v[1]**2 + v[2]**2)
	rho = math.acos(v1_dot_v2 / (mag_v1 * mag_v2))
	rho = -rho
	print(rho)
	rot2 = np.array([[1, 0, 0], [0, math.cos(rho), math.sin(rho)], [0, -math.sin(rho), math.cos(rho)]])

	# rotation 3 to z axis
	v1_dot_v2 = v[2]
	mag_v1 = 1.0
	mag_v2 = math.sqrt(v[0]**2 + v[2]**2)
	phi = math.acos(v1_dot_v2 / (mag_v1 * mag_v2))
	phi = -phi
	print(phi)
	rot3 = np.array([[math.cos(phi), 0, math.sin(phi)], [0, 1, 0], [-math.sin(phi), 0, math.cos(phi)]])


	for atom in myBGF.a:
		coord = [atom.x, atom.y, atom.z]
		newcoord = np.dot(rot1, coord)
		newcoord = np.dot(rot2, newcoord)
		newcoord = np.dot(rot3, newcoord)
		atom.x = newcoord[0]
		atom.y = newcoord[1]
		atom.z = newcoord[2]

	# save
	if isinstance(out_file, str):
		if not silent: print("Saving information to " + out_file + " ..")
		myBGF.saveBGF(out_file)
		return 1;
	else:
		return myBGF;

	##### End of pei_nnstrain


if __name__ == '__main__':

	option = ""; args = ""; bgf_file = ""; out_file = "test.bgf"; x = 0.0; y = 0.0; z = 0.0;
	usage = """
Usage: 	
"""

	if len(sys.argv) < 2:
		print(usage)
		sys.exit(0)

	# Defaults

	options, args = getopt.getopt(sys.argv[1:], 'hb:x:y:z:', ['help','bgf=','x=','y=','z='])
	for option, value in options:
		if option in ('-h', '--help'):
			print usage; sys.exit(0)
		elif option in ('-b', '--bgf'):
			bgf_file = value
		elif option in (''):
			print usage; sys.exit(0)

	#mb = bgf.BgfFile(bgf_file)
	#x1 = mb.getAtom(5)
	#x2 = mb.getAtom(17)
	#y1 = mb.getAtom(967)
	#y2 = mb.getAtom(979)
	#r1 = [(x1.x + x2.x)/2, (x1.y + x2.y)/2, (x1.z + x2.z)/2]
	#r2 = [(y1.x + y2.x)/2, (y1.y + y2.y)/2, (y1.z + y2.z)/2]
	#v = [r2[0]-r1[0], r2[1]-r1[1], r2[2]-r1[2]]
	v = [1, 0, 0]
	#print(v)
	rotate(bgf_file, v, out_file, silent=False)
