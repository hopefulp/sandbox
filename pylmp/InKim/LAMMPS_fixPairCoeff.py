#!/home/noische/Enthought/Canopy_64bit/User/bin/python
"""
LAMMPS_fixPairCoeff.py
Original: Dec 28 2011 In Kim
"""

# Python Modules
import sys
import os
import string
import time
import getopt
import math
import lammps

# Custom Modules
sys.path.append("/home/noische/scripts")
sys.path.append("/home/noische/script")
import nutils as nu

# Globals
version = '120127'

if __name__ == '__main__':

	option = ""; args = ""; in_file = ""; size = 0.0; out_file = "";
	number = 0
	usage = """
Usage: LAMMPS_fixPairCoeff.py -i lammps_in_file -o output_file

Options are:
	-i	LAMMPS input file.
	-o	Output LAMMPS input file.
	"""

	if len(sys.argv) < 2:
		print(usage); sys.exit(0)

	options, args = getopt.getopt(sys.argv[1:], 'hi:o:', ['help','in=','out='])
	for option, value in options:
		if option in ('-h', '--help'):
			print usage; sys.exit(0)
		elif option in ('-i', '--in'):
			in_file = value
		elif option in ('-o', '--out'):
			out_file = value
		elif option in (''):
			print(usage); sys.exit(0)

	# default settings
	if not out_file: out_file = os.path.basename(in_file) + "_mod"

	lammps.fixPairCoeff(in_file, out_file, silent=False)

