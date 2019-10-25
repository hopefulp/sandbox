#!/home/noische/program/python27/bin/python
"""
quaternize.py
Original: Apr 19 2011 In Kim
"""

# Python Modules
import sys
import os
import string
import random
import time
import getopt

# Custom Modules
import bgf
import nutils as nu
import bgftools
import clusterSetting

# Globals
version = '120728'

def protonate(bgf_file, out_file, log_file, silent=False):
	"""
protonate(bgf_file, out_file, log_file, silent=False):
	Protonate (performs alkylation of primary amines) the given bgf file.

	To-do:
	"""

	# initialization

	# open bgffile
	myBGF = bgf.BgfFile(bgf_file)

	# find all inserted hydrogen
	N_index = [];
	for atom in myBGF.a:
		if "H+" in atom.aName:
			N_index.append(atom.aNo)

	N_number = len(N_index);	# the number of all inserted hydrogen atoms
	if not silent: print(str(N_number) + " protonated sites are found.")

	# find the last HETATM (insertion position)
	last_hetatm_aNo = 0;
	for atom in myBGF.a:
		if atom.aTag != 1:
			last_hetatm_aNo = atom.aNo
			break;

	# add Cl atoms
	if len(N_index) > 0:
		for index, aNo in enumerate(N_index):
			if not silent: print("*** Protonated Hydrogen " + str(aNo) + " ***")

			atom_H = myBGF.getAtom(aNo);

			# new Cl atom near atom_H
			atom_Cl = bgf.BgfAtom();
			atom_Cl.x = atom_H.x + 3;
			atom_Cl.y = atom_H.y + 3;
			atom_Cl.z = atom_H.z + 3;
			atom_Cl.aName = "Cl"
			atom_Cl.rName = "ION";
			atom_Cl.chain = "A";
			atom_Cl.aTag = 1;
			atom_Cl.rNo = 0;
			atom_Cl.ffType = "Cl";
			atom_Cl.aNo = last_hetatm_aNo
			atom_Cl.charge = -1.0000;

			myBGF.addAtom(atom_Cl, last_hetatm_aNo-1)

		if not silent: print("Total charge: " + " : " + str(myBGF.charge()))

		# renumber
		myBGF.renumber()

		# save
		myBGF.REMARK.insert(0, "Ion added by " + os.path.basename(sys.argv[0]) + " by " + os.environ["USER"] + " on " + time.asctime(time.gmtime()))
		myBGF.saveBGF(out_file)
		if not silent: print("saving file " + out_file)

	else:
		if not silent: print("No amine groups to be counter-ionized"); sys.exit(0);

	if not silent: print("Done.")

	##### end of the quaternization #####


if __name__ == '__main__':

	option = ""; args = ""; bgf_file = ""; out_file = ""; log_file = "";
	usage = """
Usage: PEI_addCl.py -b bgfFile -o out_file -l log_file 
	"""

	if len(sys.argv) < 2:
		print(usage)
		sys.exit(0)

	options, args = getopt.getopt(sys.argv[1:], 'hb:o:l:', ['help','bgf=','output=','log='])
	for option, value in options:
		if option in ('-h', '--help'):
			print usage; sys.exit(0)
		elif option in ('-b', '--bgf'):
			bgf_file = value
		elif option in ('-o', '--option'):
			out_file = value
		elif option in ('-l', '--log'):
			log_file = value
		elif option in (''):
			print usage; sys.exit(0)

		if out_file == "": out_file = bgf_file[:-4] + "_ions.bgf"
		if log_file == "": log_file = out_file[:-4] + ".log"

	protonate(bgf_file, out_file, log_file, False)
