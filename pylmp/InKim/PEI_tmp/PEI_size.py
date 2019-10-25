#!/home/noische/program/python27/bin/python
"""
PEI_size.py
Original: Aug 23 2011 In Kim
"""

# Python Modules
import sys
import os
import string
import getopt
import math

# Custom Modules
import numpy
import bgf
import bgftools

# Globals
version = '110825'


def pei_size(bgf_file, silent=True):
	"""
	"""
	# Initialization
	dist = 0.0; max_dist = 0.0; aNo1 = 0; aNo2 = 0; is_string = True; n_atom = 0;

	# Open BGF
	if isinstance(bgf_file, bgf.BgfFile):
		myBGF = bgf_file
		is_string = False;
	else:
		if not silent: print("opening bgf file.. " + str(bgf_file))
		myBGF = bgf.BgfFile(bgf_file)
		is_string = True;

	### only considers HETATMs
	for atom in myBGF.a:
		if atom.aTag == 1:
			n_atom += 1;

	# for every atom
	for atom1 in myBGF.a:
		if atom1.aTag == 0:
			continue;

		if "PG" in atom1.rName:
			continue;

		if not silent:
			sys.stdout.write("\rReading atom " + str(atom1.aNo) + "/ " + str(n_atom))
			sys.stdout.flush()
		## for every atom
		for atom2 in myBGF.a:
			if atom2.aTag == 0:
				continue;

			if "PG" in atom2.rName:
				continue;
			
			### calculate distance
			dist = bgf.sqrDistance(atom1, atom2)
			### if this is the maximum
			if dist > max_dist:
				#### update maximum distance
				max_dist = dist
				#### save two atoms
				aNo1 = atom1.aNo
				aNo2 = atom2.aNo

	print("")
	if is_string:
		if not silent: print(str(bgf_file) + ": " + str(math.sqrt(max_dist)) + " A")
	else:
		if not silent: print(str(math.sqrt(max_dist)) + " A")
	#print(aNo1)
	#print(aNo2)

	##### End of pei_size


if __name__ == '__main__':

	option = ""; args = ""; bgf_file = ""; ff_file = ""; probability = 0; out_file = "";
	usage = """
Usage: 	
"""

	if len(sys.argv) < 2:
		print(usage)
		sys.exit(0)

	# Defaults
	safemode = False; silent = True;

	options, args = getopt.getopt(sys.argv[1:], 'hb:', ['help','bgf='])
	for option, value in options:
		if option in ('-h', '--help'):
			print usage; sys.exit(0)
		elif option in ('-b', '--bgf'):
			bgf_file = value
		elif option in (''):
			print usage; sys.exit(0)

	pei_size(bgf_file, False)
