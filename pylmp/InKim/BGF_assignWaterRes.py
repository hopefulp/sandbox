#!/opt/applic/epd/bin/python

# system modules
import sys
import os
import string
import getopt
import random

# bgf modules
sys.path.append("/home/noische/script")
import nutils as nu
import bgf
import bgftools

option = ""; args = ""; bgf_file = ""; trj_file = ""; out_file = ""; timestep = 0;
usage = """
renumberResNum.py: renumbers the residue number in the BGF file.

Usage: renumberResNum.py -b bgf_file -o out_file
"""
version = "150623"

def assignWaterResidue(bgf_file, out_file, silent=True):

	mybgf = bgf.BgfFile(bgf_file)

	for atom in mybgf.a:
		for i in bgftools.is_water(mybgf, atom.aNo):
			mybgf.getAtom(i).rName = "WAT"

	mybgf.saveBGF(out_file)

	### end of assignWaterResidue


if __name__ == '__main__':

	if len(sys.argv) < 2:
		print(usage);
		sys.exit(0)

	options, args = getopt.getopt(sys.argv[1:], 'hb:t:n:o:', ['help','bgf=','trj=','number=','out='])
	for option, value in options:
		if option in ('-h', '--help'):
			print(usage)
			sys.exit(0);
		elif option in ('-b', '--bgf'):
			bgf_file = value
		elif option in ('-o', '--out'):
			out_file = value
		elif option == NULL:
			print(usage)
			sys.exit(0)

	# default
	if out_file == "":
		out_file = os.path.basename(bgf_file).split(".bgf")[0] + "_resname.bgf"


	# main call
	print(sys.argv[0] + " version " + str(version))
	assignWaterResidue(bgf_file, out_file, silent=False)

