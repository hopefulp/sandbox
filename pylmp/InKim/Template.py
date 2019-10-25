#!/home/noische/program/python27/bin/python
"""
template.py
Original: Dec 28 2011 In Kim
"""

# Python Modules
import sys
import os
import string
import random
import time
import getopt

# Custom Modules
sys.path.append("/home/noische/scripts")
sys.path.append("/home/noische/script")
import bgf
import bgftools
import nutils as nu

# Globals
version = '111228'

def template(bgf_file, out_file, silent=False):
	"""
def template():
	Write what this function does.

Function Parameters:
	bgf_file	A string of filename or BgfFile class.
	out_file	A string of filename or BgfFile class.
	"""
	# open BGF
	if isinstance(bgf_file, bgf.BgfFile):
		myBGF = bgf_file
	else:
		if not silent: print("opening bgf file.. " + str(bgf_file))
		myBGF = bgf.BgfFile(bgf_file)

	# Write down the Template function here


	# save BGF
	if isinstance(out_file, str):
		if not silent: print("Saving information to " + out_file + " ..")
		myBGF.saveBGF(out_file)
		return 1;
	else:
		return myBGF;


	### end of template


if __name__ == '__main__':

	option = ""; args = ""; bgf_file = ""; size = 0.0; out_file = "";
	number = 0
	usage = """
Usage: template.py -b "bgfFiles" -s size -n monomers -o output

Options are:
	-b	Input BGF file.
	-o	Output BGF file.
	"""

	if len(sys.argv) < 2:
		print(usage); sys.exit(0)

	options, args = getopt.getopt(sys.argv[1:], 'hb:o:', ['help','bgf=','output='])
	for option, value in options:
		if option in ('-h', '--help'):
			print usage; sys.exit(0)
		elif option in ('-b', '--bgf'):
			bgf_file = value
		elif option in ('-o', '--output'):
			out_file = value
		elif option in (''):
			print(usage); sys.exit(0)

	# default settings
	if not out_file: out_file = os.path.basename(bgf_file).split(".bgf")[0] + "_mod" + ".bgf"

	template(bgf_file, out_file, silent=False)

