#!/home/noische/python

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

option = ""; args = ""; bgf_file = ""; trj_file = ""; out_file = ""; timestep = 0;
usage = """
countMoleculesBGF.py: counts the number of molecules in the BGF file.

Usage: countMoleculesBGF.py -b bgf_file 
"""
version = "111007"

"""
* Updates
- 111007: Counts the number of water molecules
"""

def countMoleculeNum(bgf_file, silent=True):

	n_water = 0;

	l_molecule = bgftools.getMoleculeList(bgf_file)
	natom = len(nu.flatten(l_molecule))
	nmol = len(l_molecule)
	l_molecule_atoms = [];
	for cluster in l_molecule:
		l_molecule_atoms.append(len(cluster))
		if len(cluster) == 3:
			n_water += 1

	if not silent: print(str(nmol) + " Molecules (" + str(natom) + " atoms) exists in the BGF file.")
	if not silent: print("Number of water molecules (i.e. natoms = 3): " + str(n_water))
	#if not silent: print("Size of molecules: " + str(l_molecule_atoms))

	return nmol

	### end of countMoleculeNum


if __name__ == '__main__':

	if len(sys.argv) < 2:
		print(usage);
		sys.exit(0)

	options, args = getopt.getopt(sys.argv[1:], 'hb:', ['help','bgf='])
	for option, value in options:
		if option in ('-h', '--help'):
			print(usage)
			sys.exit(0);
		elif option in ('-b', '--bgf'):
			bgf_file = value
		elif option == NULL:
			print(usage)
			sys.exit(0)

	# main call
	print(sys.argv[0] + " version " + str(version))
	countMoleculeNum(bgf_file, silent=False)

