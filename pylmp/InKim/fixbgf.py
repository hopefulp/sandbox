#!/home/noische/program/python27/bin/python
"""
fixbgf.py
Original: Apr 20 2011 In Kim
"""

# Python Modules
import sys
import os
import string
import getopt
import time

# Custom Modules
import bgf

# Globals
version = '110420'


def fixbgf(bgf_file, out_file, silent=False):
	"""
fixbgf(bgf_file, out_file):

	bgf_file	A string for input file.
	out_file	A string for output file.

	To-do:
		None
	"""
	# Open BGF
	myBGF = bgf.BgfFile(bgf_file)

	# Save BGF
	myBGF.REMARK.insert(0, "Renumbered by " + os.path.basename(sys.argv[0]) + " by " + os.environ["USER"] + " on " + time.asctime(time.gmtime()))
	myBGF.renumber()
	myBGF.saveBGF(out_file)

	##### End of fixbgf


if __name__ == '__main__':

	option = ""; args = ""; bgf_file = ""; ff_file = ""; probability = 0; out_file = "";
	usage = """
Usage: fixbgf.py -b bgfFile -o outFile 
	This script fixes small bugs in a Biograf BGF file generated from CERIUS2.

Options are:
	-b	REQUIRED. An input BGF file
	-o	OPTIONAL. An output BGF file. The suffix "_fixed" will be attached if not stated.

	Report any bugs to in.kim@kaist.ac.kr
	"""

	if len(sys.argv) < 2:
		print(usage)
		sys.exit(0)

	options, args = getopt.getopt(sys.argv[1:], 'hb:o:', ['help','bgf=','output='])
	for option, value in options:
		if option in ('-h', '--help'):
			print usage; sys.exit(0)
		elif option in ('-b', '--bgf'):
			bgf_file = value
		elif option in ('-o', '--option'):
			out_file = value
		elif option in (''):
			print usage; sys.exit(0)

		if out_file == "": out_file = bgf_file[:-4] + "_fixed.bgf"

	fixbgf(bgf_file, out_file, silent=False)
