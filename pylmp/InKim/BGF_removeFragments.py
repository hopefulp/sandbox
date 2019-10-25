#!/home/noische/program/python27/bin/python

# system modules
import sys
import os
import string
import getopt
import random

# bgf modules
sys.path.append("/home/noische/script")
import nutils as nu
import bgftools

option = ""; args = ""; bgf_file = ""; n_frag = 0; out_file = ""; resname = "";
usage = """
removeFragments.py: deletes molecules in BGF which have less atoms than requested.

Usage: removeFragments.py -b bgf_file -n atoms -o out_file
"""
version = "110817"

def removeFragments(bgf_file, resname, n_frag, out_file, silent=True):

	bgftools.removeFragments(bgf_file, resname, n_frag, out_file, silent=False)

	### end of countMoleculeNum


if __name__ == '__main__':

	if len(sys.argv) < 2:
		print(usage);
		sys.exit(0)

	options, args = getopt.getopt(sys.argv[1:], 'hb:n:o:r:', ['help','bgf=','size=','output=','residue='])
	for option, value in options:
		if option in ('-h', '--help'):
			print(usage)
			sys.exit(0);
		elif option in ('-b', '--bgf'):
			bgf_file = value
		elif option in ('-n', '--size'):
			n_frag = int(value)
		elif option in ('-o', '--output'):
			out_file = value
		elif option in ('-r', '--residue'):
			resname = value
		elif option == NULL:
			print(usage)
			sys.exit(0)
		else:
			print(usage)
			sys.exit(1)

	# main call
	print(sys.argv[0] + " version " + str(version))
	removeFragments(bgf_file, resname, n_frag, out_file, silent=False)

