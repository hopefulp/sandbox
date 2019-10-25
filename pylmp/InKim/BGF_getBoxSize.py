#!/opt/applic/epd/bin/python

import sys
import re
import string
import getopt
import optparse
import math
import time

from os import popen

#-----------------
# update the coordinate in the original BGF file from LAMMPS trajectory file
#_________________
def getBgfBoxSize(bgf_file):
	boxsize = [0, 0, 0, 0, 0, 0]

	f_bgf_file = open(bgf_file)

	# parse the original bgf and write to outfile
	while 1:
		line = f_bgf_file.readline()
		if not line:
			break
		if 'CONECT' in line:
			break

		# very simple parsing method.. just split and get it back again with just updating
		if 'HETATM' in line:
			line = f_bgf_file.readline()
			parse = re.split('\s*', line)
			x = float(parse[6]) # scaled x coord
			y = float(parse[7]) # scaled y coord
			z = float(parse[8]) # scaled z coord

			if x < boxsize[0]:
				boxsize[0] = x
			elif x > boxsize[1]:
				boxsize[1] = x

			if y < boxsize[2]:
				boxsize[2] = y
			elif y > boxsize[3]:
				boxsize[3] = y

			if z < boxsize[4]:
				boxsize[4] = z
			elif z > boxsize[5]:
				boxsize[5] = z

	a = []
	a.append(boxsize[1]-boxsize[0])
	a.append(boxsize[3]-boxsize[2])
	a.append(boxsize[5]-boxsize[4])
	#print(max(a))
	print("Range of " + bgf_file + ": " + str(boxsize[0]) + " < x < " + str(boxsize[1]) + ", " + str(boxsize[2]) + " < y < " + str(boxsize[3]) + ", " + str(boxsize[4]) + " < z < " + str(boxsize[5]))
	print("Size of  " + bgf_file + ": x = " + str(boxsize[1]-boxsize[0]) + " y = " + str(boxsize[3]-boxsize[2]) + " z = " + str(boxsize[5]-boxsize[4]))

	return 1

if __name__ == "__main__":

	option = ""; args = ""; bgf_file = ""
	usage = """
Usage: python getBgfBoxSize.py -b bgf_file
	Note that this script does not consider vdw radius of atoms.
	"""

        if len(sys.argv) < 2:
                print(usage)
		sys.exit(0)

	options, args = getopt.getopt(sys.argv[1:], 'hb:t:o:', ['help','bgf='])
	for option, value in options:
	        if option in ('-h', '--help'):
	                print(usage); sys.exit(0)
	        elif option in ('-b', '--bgf'):
	                bgf_file = value
	        elif option in NULL:
	                print(usage); sys.exit(0)

	# main call
	getBgfBoxSize(bgf_file)
