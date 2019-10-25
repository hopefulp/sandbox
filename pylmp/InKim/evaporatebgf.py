#!/home/noische/program/python27/bin/python
"""
evaporatebgf.py
Original: Apr 27 2011 In Kim
"""

# Python Modules
import sys
import os
import string
import getopt
import time

# Custom Modules
import bgf
import bgftools

# Globals
version = '110427'


def evaporate(bgf_file, out_file=False, safemode=False, silent=True):
	"""
evaporate(bgf_file, out_file):

	bgf_file	A string for input file.
	out_file	A string for output file.
	safemode	A boolean for safe mode determination.

	To-do:
		None
	"""
	# Initialization
	temp = [];

	# Open BGF
	if isinstance(bgf_file, bgf.BgfFile):
		myBGF = bgf_file
	else:
		if not silent: print("opening bgf file.. " + str(bgf_file))
		myBGF = bgf.BgfFile(bgf_file)

	# Delete water molecules
	if safemode:
		for atom in myBGF.a:
			temp.append(atom.aNo)
		for aNo in temp:
			bgftools.deleteWaterAtoms(myBGF, aNo)
	else:
		if not silent: print("listing water molecules..")
		for atom in myBGF.a:
			if "WAT" in atom.rName:
				temp.append(myBGF.a2i[atom.aNo])
		n_delatoms = len(temp)
		if not silent: print("found " + str(n_delatoms) + " atoms for water molecules..")
                if n_delatoms % 3 != 0:
			nu.die("ERROR: the number of atoms to be deleted is not divided in 3. Check the structure.")
		if not silent: print("deleting water molecules..")
		myBGF.delAtoms(temp, silent)
		myBGF.renumber()

	# Save BGF
	if isinstance(out_file, str):
		if not silent: print("\nsaving information to " + out_file + " ..")
		myBGF.REMARK.insert(0, "Evaporated by " + os.path.basename(sys.argv[0]) + " by " + os.environ["USER"] + " on " + time.asctime(time.gmtime()))
		myBGF.saveBGF(out_file)
		return 1;
	else:
		return myBGF;

	##### End of evaporate


if __name__ == '__main__':

	option = ""; args = ""; bgf_file = ""; ff_file = ""; probability = 0; out_file = "";
	usage = """
Usage: evaporatebgf.py -b bgfFile -o outFile 
	This script evaporates all water molecules in a BGF file.
	If a residue name is same as "WAT", then the atom will be deleted.

Options are:
	-b	REQUIRED. An input BGF file
	-o	OPTIONAL. An output BGF file. The suffix "_fixed" will be attached if not stated.
	-s	OPTIONAL. Activates safe mode. 
		By default, all atoms which have the residue name "WAT" will be deleted.
		In a safe mode, the script will automatically check whether it is a water molecule or not.
		

	Report any bugs to in.kim@kaist.ac.kr
	"""

	if len(sys.argv) < 2:
		print(usage)
		sys.exit(0)

	# Defaults
	safemode = False; silent = True;

	options, args = getopt.getopt(sys.argv[1:], 'hb:o:sv', ['help','bgf=','output=','safe=','verbose='])
	for option, value in options:
		if option in ('-h', '--help'):
			print usage; sys.exit(0)
		elif option in ('-b', '--bgf'):
			bgf_file = value
		elif option in ('-o', '--option'):
			out_file = value
		elif option in ('-s', '--safe'):
			safemode = True
		elif option in ('-v', '--verbose'):
			silent = False
		elif option in (''):
			print usage; sys.exit(0)

		if out_file == "": out_file = bgf_file[:-4] + "_dry.bgf"

	evaporate(bgf_file, out_file, safemode, silent=False)
