#!/opt/applic/epd/bin/python

import sys
import re
import string
import getopt
import optparse
import math
import os
import time

import bgf
import nutils as nu

def bgf2qcin(bgf_file, qchem_in_file, rem_file):

	# bgf open
	print("Reading " + bgf_file + " ..")
	myBGF = bgf.BgfFile(bgf_file)


	# output
	output = ["$comment\n", "from " + bgf_file + " at " + time.asctime(time.gmtime()) + " on " + os.environ["HOSTNAME"] + " by " + os.environ["USER"] + "\n", "$end\n"]

	# charge and multiplicity
	output += ["\n$molecule\n"]
	chg = myBGF.charge()
	if abs(chg) < 0.000001:
		output += ["0" + "\t" + "1" + "\n"]
	else:
		nu.warn("Charge is not neutral. Charge will be written as " + "{0:<4.1f}".format(chg))
		output += ["{0:<4.1f}".format(chg) + "\t" + "1" + "\n"]

	# coordinates
	for atom in myBGF.a:
		output += [atom.ffType.split("_")[0] + "\t" + "{0:10.6f} {1:10.6f} {2:10.6f}".format(atom.x, atom.y, atom.z) + "\n"]

	output += ["$end\n\n"]

	# rem
	remf = open(rem_file)
	remline = remf.readlines()
	output += remline

	# target file
	qcinf = open(qchem_in_file, 'w')
	qcinf.writelines(output)
	qcinf.close()

	print("Generated " + qchem_in_file + " . Done.")


	### end of function


if __name__ == "__main__":
	bgf_file = ""; qchem_in_file = ""; rem_file = ""; 
	usage = """
BGF_Bgf2QchemInput.py -b bgf_file -o qchem_in_file -r rem_file
	"""


	options, args = getopt.getopt(sys.argv[1:], 'hb:o:r:', ['help', 'bgf=', 'qcin=', 'rem='])

	if len(sys.argv) < 2:
		print(usage)
		sys.exit(0)

	print("Requested options: " + str(options))

	for option, value in options:
	        if option in ('-h', '--help'):
	                print(usage)
			sys.exit(0)
	        elif option in ('-b', '--bgf'):
	                bgf_file = value
	        elif option in ('-o', '--qcin'):
	                qchem_in_file = value
	        elif option in ('-r', '--rem'):
	                rem_file = value
	        elif option == NULL:
			print(usage)
			sys.exit(0)

	# required options
	if bgf_file == "" or os.path.exists(bgf_file) == False:
		nu.die("No BGF file named " + bgf_file)

	# default options
	if qchem_in_file == "":
		qchem_in_file = bgf_file.split(".bgf")[0] + ".in"
	
	# main call
	bgf2qcin(bgf_file, qchem_in_file, rem_file)
