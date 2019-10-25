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
R = 0.0;	# radius
sOC = 0.0; 	# sigma_OC parameter
L = 0.0; 	# CNT length
pi = math.pi;	# pi
u = 0.0; 
prefix = "";

usage = """
countWaterCNT.py -b bgf_file -t trj_file -n #step
"""
version = "130729"

#-----------------
# 
# 
# 
#_________________

def sortkey(list):
	return list[0];

def getSlipLength(force_profile, vel_profile, silent=False):
	"""
	This calculates the profile of averaged velocity according to r. 
	NOTE: This is different from averaged velocity profile according to r.
	"""

	### init
	line = []; n_header = 0; output = "";
	output += "t\tforce\tarea\tvz\tvisc\tsliplen\n"

	myOUT = open(str(prefix) + ".sliplength.profile", 'w')


	### read force data from force profile
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


	### read vel file
	if not silent: print("Reading velocity profile..")
	f_vp = open(vel_profile)
	vel = dict();
	while 1:
		line = f_vp.readline()
		if not line:
			break;
		if line[0] == "#":
			continue;
		parse = line.split()
		step = int(parse[0])
		vel[step] = float(parse[1])

	k = force.keys();
	k.sort()


	### calculate slip length
	for i in k:

		v_slip = vel[i] * 1e-10 / 1e-15	# A/fs -> m/s
		f = force[i] * (6.9478e-21) / (1e-10)	# kcal/mol.A -> J/m
		reff = R - 0.5 * sOC
		A = 2 * pi * reff * L
		A *= (1e-20)		# A2 -> m2
		l = f / A / v_slip	# lambda
		b = u / l * 1e9		# slip length in nm

		output += str(i) + "\t" + str(force[i]) + "\t" + str(A) + "\t" + str(v_slip) + "\t" + str(u) + "\t" + str(l) + "\t" + str(b) + "\n"


	myOUT.write(output)
	myOUT.close()

	print('')
	return 1

	### end of function


if __name__ == "__main__":
	trj_file = ""; force_profile = ""; vel_profile = "";
	options, args = getopt.getopt(sys.argv[1:], 'hi:f:v:R:r:L:u:s:', ['help', 'force=', 'vel=', 'R=', 'r=', 'L=', 'viscosity=', 'prefix='])

	if len(sys.argv) < 2:
		print(usage)
		sys.exit(0)

	print "Requested options: " + str(options)

	for option, value in options:
	        if option in ('-h', '--help'):
	                print(usage)
			sys.exit(0)
		elif option in ('-f', '--force'):
			force_profile = value
		elif option in ('-v', '--vel'):
			vel_profile = value
	        elif option in ('-R', '--R'):
			R = float(value)
		elif option in ('-r', '--sigma_OC='):
			sOC = float(value)
		elif option in ('-L', '--length'):
			L = float(value)
		elif option in ('-u', '--viscosity'):
			u = float(value);	# in kg/m.s unit
		elif option in ('-s', '--prefix'):
			prefix = value
	        elif option == NULL:
			print(usage)
			sys.exit(0)
	
	# more options: resname CNT, resname solvent, save BGF

	# main call
	getSlipLength(force_profile, vel_profile)
