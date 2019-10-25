#!/home/noische/python

# Python Modules
import sys
import re
import string
import getopt
import optparse
from os import popen
from types import *

# BGF modules
sys.path.append("/home/noische/script")
import bgf
import bgftools
import nutils as nu

# Globals
option = ""; args = ""; bgf_file = ""; thresh = "";
usage = """
removeBadContacts.py: read coordinate data from the LAMMPS trajectory file
           write the data to the original BGF file
Usage: removeBadContacts.py -b bgf_file -t thresh -o out_file
"""

#-----------------
# update the coordinate in the original BGF file from LAMMPS trajectory file
#_________________
def removebadcontacts(bgf_file, out_file, thresh, silent=True):

	if isinstance(bgf_file, bgf.BgfFile):
		myBGF = bgf_file
	else:
		if not silent: print(sys.argv[0] + ": Removing bad contacts from " + bgf_file + " with distance threshold " + str(thresh) + " A and saving to " + out_file + ".")
		myBGF = bgf.BgfFile(bgf_file);

	delete_list = [];

	if not silent: print("removeBadContacts will remove water molecules within " + str(thresh) + " Angstrom from the solute." )
	for solute_aNo in bgftools.listSoluteAtoms(myBGF):
		solute = myBGF.getAtom(solute_aNo)
		for oxygen_aNo in bgftools.listOxygenAtoms(myBGF):
			water = bgftools.is_water(myBGF, oxygen_aNo)
			if water != [] or type(water) != NoneType:
				if len(water) == 3:		# if water is a water molecule,
					# calculate the distance between solute atoms and the molecule
					O = myBGF.getAtom(water[0])
					H1 = myBGF.getAtom(water[1])
					H2 = myBGF.getAtom(water[2])
					dist_O_sol = bgf.distance(solute, O)
					dist_H1_sol = bgf.distance(solute, H1)
					dist_H2_sol = bgf.distance(solute, H2)
					if dist_O_sol < thresh or dist_H1_sol < thresh or dist_H2_sol < thresh:
						if oxygen_aNo not in delete_list:
							delete_list.append(oxygen_aNo)
							delete_list.sort()
						##if water not in delete_list:
							##delete_list.append(water)

	if delete_list != []:
		delete_list.reverse()	# reverse sort for delAtoms
		for oxygen_index in delete_list:
			bgftools.deleteWaterAtoms(myBGF, oxygen_index)
		myBGF.renumber()

		if not silent: print("removeBadContacts: " + str(len(delete_list)) + " water molecules are removed.")
	else:
		if not silent: print("There are no water molecules that corresponds to the criteria.")

	if isinstance(out_file, bgf.BgfFile):
		return myBGF;
	else:
		myBGF.saveBGF(out_file)
		return 1;

	### end of removeBadContacts.py

if __name__ == "__main__":

	options, args = getopt.getopt(sys.argv[1:], 'hb:t:o:', ['help','bgf=','thresh=','out='])
	print("Requested options: " + str(options))
	for option, value in options:
		if option in ('-h', '--help'):
			print(usage)
			sys.exit(0)
		elif option in ('-b', '--bgf'):
			bgf_file = value
		elif option in ('-t', '--thresh'):
			thresh = float(value)
		elif option in ('-o', '--out'):
			out_file = value
		elif option == NULL:
			print(usage)
			sys.exit(0)

	# main call
	removebadcontacts(bgf_file, out_file, thresh, False)


