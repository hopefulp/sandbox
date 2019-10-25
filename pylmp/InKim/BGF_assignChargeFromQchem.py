#!/opt/applic/epd/bin/python

import sys
import re
import string
import getopt
import optparse
import math
import time

import bgf
import nutils as nu

version = '130111'

def listRightIndex(alist, value):
	"""
	Finds a biggest index of the given value in the list.
	Reference: http://stackoverflow.com/questions/9836425/equivelant-to-rindex-for-lists-in-python (EwyynTomato)
	"""
	return len(alist) - alist[-1::-1].index(value) - 1


#-----------------
# read charges from qchem output and assign to bgf
#_________________
def assignCharges(bgf_file, qcout_file, out_file, silent=True):

	# read qcout file
	f_qcout_file = open(qcout_file)
	qcdata = [];

	while 1:
		line = f_qcout_file.readline();

		if not line:
			break;

		qcdata.append(line.strip())

	# find 'MISSION COMPLETE'
	if not silent: print("Checking whether Q-Chem calculation is properly completed.")
	fileTester = False;
	for i in qcdata:
		if "MISSION COMPLETED" in i:
			fileTester = True;
			break;
	if fileTester != True:
		nu.warn("Q-Chem output file is incomplete.")

	# find the last "Ground-State Mulliken Net Atomic Charges"
	i1 = listRightIndex(qcdata, "Ground-State Mulliken Net Atomic Charges")
	qcdata = qcdata[i1:]

	# if not, die
	if i1 == 0:
		nu.die("No Mulliken Charge Information in " + qcout_file)

	# find "Sum of atomic charges"
	for index, i in enumerate(qcdata):
		if "Sum of atomic charges" in i:
			i2 = index;

	qcdata = qcdata[:i2]

	# get Charges part
	temp = [];
	i1 = qcdata.index('----------------------------------------')
	i2 = listRightIndex(qcdata, '----------------------------------------')
	qcdata = qcdata[i1 + 1: i2]
	for i in qcdata:
		parse = re.split('\s*', i)
		temp.append(parse)
	qcdata = temp;
	if not silent: print("Found the Charge Information.")

	# open bgf
	if isinstance(bgf_file, bgf.BgfFile):
		myBGF = bgf_file
	else:
		if not silent: print("Reading " + bgf_file + " ..")
		myBGF = bgf.BgfFile(bgf_file)

	# check the number of atoms in bgf
	n_atoms_in_qcout = len(qcdata)
	n_atoms_in_bgf = len(myBGF.a)
	if n_atoms_in_qcout != n_atoms_in_bgf:
		nu.die("Number of atoms in both files are different. \n\t" + qcout_file + ": " + str(n_atoms_in_qcout) + " atoms and " + bgf_file + ": " + str(n_atoms_in_bgf) + " atoms.")

	# check the order of atoms in both files and assign charges if right
	for index, i in enumerate(qcdata):
		atom = myBGF.a[index]
		if i[1] in atom.ffType:
			pass;
		else:
			nu.die("The order of the atoms are different: Check if the Q-Chem input file is directly generated from the BGF file.")

	# update charges
	if not silent: print("Updating Charges..")
	for index, i in enumerate(qcdata):
		atom = myBGF.a[index]
		atom.charge = float(i[2])

	# save bgf
	if isinstance(out_file, str):
		if not silent: print("Saving information to " + out_file + " ..")
		myBGF.saveBGF(out_file)
		return 1;
	else:
		return myBGF;


if __name__ == "__main__":

	option = ""; args = ""; bgf_file = ""; qcout_file = ""; out_file = "";
	usage = """
Usage: BGF_assignChargeFromQchem.py -b bgf_file -q qcout_file -o out_file
	"""

        if len(sys.argv) < 2:
                print(usage)
		sys.exit(0)

	options, args = getopt.getopt(sys.argv[1:], 'hb:q:o:', ['help','bgf=','qcout=','out='])
	for option, value in options:
	        if option in ('-h', '--help'):
	                print(usage); sys.exit(0)
	        elif option in ('-b', '--bgf'):
	                bgf_file = value
		elif option in ('-q', '--qcout'):
			qcout_file = value
		elif option in ('-o', '--out'):
			out_file = value
	        elif option in NULL:
	                print(usage); sys.exit(0)

	# main call
	print(sys.argv[0] + " version " + str(version))
	assignCharges(bgf_file, qcout_file, out_file, False)
