#!/home/joonho/anaconda3/bin/python

"""
addSolvent.py
Original: Mar 27 2011 In Kim
Version: 110514

# Version updates:
- 110514: Solvent box is extended to general solvent cube 'solvent_file' from f3c_water box equilibrated by Tod Pascal.

# TODO: Extend remove bad contacts with other solvent molecules.
"""

# Python Modules
import sys
import os
import string
import random
import time
import math

# Custom Modules
sys.path.append("/home/noische/script")
sys.path.append("/home/noische/scripts")
import bgf
import nutils as nu
from dump import *
import bgftools
from removeBadContacts import *

# Globals
version = '120924'


def addsolvent(bgf_file, solvent_bgf, size, margin, out_file, ff_file, silent=True):

	### initialize
	water = False;

	### load the solute bgf file
	if not silent: print("Initializing..")
	myBGF = bgf.BgfFile(bgf_file)	# str

	#solventBGF = bgf.BgfFile("/home/noische/scripts/dat/WAT/f3c_box.bgf")	# F3C waterbox
	if not silent: print("Loading the solvent file " + solvent_bgf + " ..")
	solventBGF = bgf.BgfFile(solvent_bgf)

	### Generate error when the solvent box is not periodic:
	if not solventBGF.PERIOD:
		nu.die("addSolvent: The solvent file is not periodic. Use a box full of solvent.")

	### Check the type of solvent
	if not silent: print("(the solvent box seems to be full of " + os.path.basename(solvent_file)[:-4] + " )")
	if "f3c" in solvent_bgf: water = True		# this flag is used to remove the bad contacts with the molecule.
	if "spc" in solvent_bgf: water = True		# this flag is used to remove the bad contacts with the molecule.
	if "tip" in solvent_bgf: water = True		# this flag is used to remove the bad contacts with the molecule.

	### calculate the box size
	if not silent: print("Analyzing box information..")
	if len(myBGF.CRYSTX) > 2:
		strsize = [ myBGF.CRYSTX[0], myBGF.CRYSTX[1], myBGF.CRYSTX[2] ]
	else:
		box = bgf.getBGFSize(myBGF, 0)          # [xlo, xhi, ylo, yhi, zlo, zhi]
		strsize = [ box[1] - box[0], box[3] - box[2], box[5] - box[4] ]

	###
	if size == "" and margin != "":
		strsize = [ strsize[0] + 2 * margin[0], strsize[1] + 2 * margin[1], strsize[2] + 2 * margin[2] ]	# add margin on str size
	elif size != "" and margin == "":
		strsize = size

	waterboxsize = solventBGF.CRYSTX[:3]	# REMARK: This is a kind of constant.
	copyNumber = [0, 0, 0];
	for index, i in enumerate(copyNumber):
		copyNumber[index] = math.ceil(strsize[index] / waterboxsize[index])	# how many times to replicate the water box

	if not silent: print("Creating box information: " + str(strsize))

	### replicate the solvent box
	bigboxBGF = bgftools.replicateCell(solventBGF, copyNumber, True)
	bigboxBGF.saveBGF("_replicate.bgf")
	if not silent: print("- Number of atoms in the created box: " + str(len(bigboxBGF.a)))

	### trim the water box
	if water:
		if not silent: print("Generating water box.. Calculating water molecules")
		delatom = []; delwater = []; delwaterindex = [];
		for atom in bigboxBGF.a:
			if atom.x > strsize[0] or atom.y > strsize[1] or atom.z > strsize[2]:
				delatom.append(atom.aNo)
		for aNo in delatom:
			water_molecule = bgftools.is_water(bigboxBGF, aNo)
			if not water_molecule in delwater:
				delwater.append(water_molecule)
		delwater = nu.flatten(delwater)
		delwater = nu.removeRepeat(delwater)
		delwater.sort()
		delwater.reverse()
		for aNo in delwater:
			delwaterindex.append(bigboxBGF.getAtomIndex(aNo))
		del(delwater)
		if not silent: print("Generating water box.. Trimming")
		bigboxBGF.delAtoms(delwaterindex, False)
		bigboxBGF.renumber()
	elif not water:
		if not silent: print("Generating solvent box.. Extracting solvent molecules")
		delatom = []; delsolvent = []; delsolventindex = [];
		for atom in bigboxBGF.a:
			if atom.x > strsize[0] or atom.y > strsize[1] or atom.z > strsize[2]:
				delatom.append(atom.aNo)
		for aNo in delatom:
			molecule_list = []
			molecule = bgftools.getmolecule(bigboxBGF, bigboxBGF.getAtom(aNo), molecule_list)
			for number in molecule_list:
				if not number in delsolvent: delsolvent.append(number)
		delsolvent = nu.flatten(delsolvent)
		delsolvent = nu.removeRepeat(delsolvent)
		delsolvent.sort()
		delsolvent.reverse()
		for aNo in delsolvent:
			delsolventindex.append(bigboxBGF.getAtomIndex(aNo))
		if not silent: print("Generating solvent box.. Trimming")
		bigboxBGF.delAtoms(delsolventindex, False)
		bigboxBGF.renumber()

	### merge two structure
	# REMARK: it is natural to have the periodic information of water box for the output BGF file.
	# REMARK: HETATOM should be located on the first of the BGF file. So use dummy for merging.
	if not silent: print("\nAdding trimmed solvent box to the structure..")
	bigboxcenter = [ strsize[0] / 2, strsize[1] / 2, strsize[2] / 2]
	a, b, c = bgftools.getCom(myBGF, ff_file)
	bgf.moveBGF(myBGF, bigboxcenter[0] - a, bigboxcenter[1] - b, bigboxcenter[2] - c)

	## remove bad contacts between solutes and solvents
	if not silent: print("Atom distance calculation for contacts..")
	delatom = []; delsolvent = []; delsolventindex = [];
	for atom1 in myBGF.a:
		for atom2 in bigboxBGF.a:
			# if the distance between atom1 and atom2 is less than 2.8, add to a delete list
			dist_sq = (atom1.x - atom2.x)**2 + (atom1.y - atom2.y)**2 + (atom1.z - atom2.z)**2
			if dist_sq < 7.84:
				delatom.append(atom2.aNo)

	# delete bad atoms!
	if not silent: print("Removing bad contacts..")
	for aNo in delatom:
		molecule_list = [];
		molecule = bgftools.getmolecule(bigboxBGF, bigboxBGF.getAtom(aNo), molecule_list)
		for number in molecule_list:
			if not number in delsolvent: delsolvent.append(number)
	delsolvent = nu.flatten(delsolvent)
	delsolvent = nu.removeRepeat(delsolvent)
	delsolvent.sort(); delsolvent.reverse();
	for aNo in delsolvent:
		delsolventindex.append(bigboxBGF.getAtomIndex(aNo))
	bigboxBGF.delAtoms(delsolventindex, False)
	bigboxBGF.renumber()

	## compute stats for adding solvents
	if not silent: print("\nComputing stats..")
	mol_list = bgftools.getMoleculeList(bigboxBGF)
	n_mol = len(mol_list)
	n_atom = len(nu.flatten(mol_list))
	if not silent: print(str(n_mol) + " molecules (" + str(n_atom) + " atoms) will be added.")

	## merge
	bigboxBGF = myBGF.merge(bigboxBGF, True)
	if not silent: print("Total atoms in the file: " + str(len(bigboxBGF.a)))

	## some paperworking for periodic box
	bigboxBGF.OTHER = solventBGF.OTHER
	bigboxBGF.PERIOD = solventBGF.PERIOD
	bigboxBGF.AXES = solventBGF.AXES
	bigboxBGF.SGNAME = solventBGF.SGNAME
	bigboxBGF.CRYSTX = solventBGF.CRYSTX
	bigboxBGF.CELLS = solventBGF.CELLS

	## adjust the size of the box
	bigboxBGF.CRYSTX = [ strsize[0], strsize[1], strsize[2], solventBGF.CRYSTX[3], solventBGF.CRYSTX[4], solventBGF.CRYSTX[5] ]

	"""
	### remove bad contacts and save
	if water:
		if not silent: print("Removing bad contacts: distance criteria is 2.8 A")
		bigboxBGF = removebadcontacts(bigboxBGF, bigboxBGF, 2.8)

		### renumber residue numbers
		# RULE: all rNos for water molecules will be renumbered from 2. (for createLammpsInput.pl)
		if not silent: print("Renumbering water molecules..")
		max_rNo_in_hetatm = 0;
		oxygen_list = bgftools.listOxygenAtoms(bigboxBGF)
	
		### find the biggest rNo among HETATM
		for atom in bigboxBGF.a:
			if atom.aTag == 1:
				if max_rNo_in_hetatm < atom.rNo: max_rNo_in_hetatm = atom.rNo
	
		### update rNo in water molecules
		rNo_for_water = max_rNo_in_hetatm + 500;
		for aNo in oxygen_list:
			water_aNo = bgftools.is_water(bigboxBGF, aNo)	# get aNo in a water molecule by checking the oxygen atom
			if water_aNo != []:
				for atom_aNo in water_aNo:
					bigboxBGF.getAtom(atom_aNo).rNo = rNo_for_water
				rNo_for_water += 1
	"""
	
	### record BGF remarks
	bigboxBGF.REMARK.insert(0, "Solvents added by " + os.path.basename(sys.argv[0]) + " by " + os.environ["USER"] + " on " + time.asctime(time.gmtime()))
	bigboxBGF.REMARK.insert(0, "Solvents: " + str(solvent_file))

	### save BGF
	if not silent: print("Saving the file.. see " + str(out_file))
	bigboxBGF.saveBGF(out_file)

	return 1;

	### end of addSolvent.py


if __name__ == '__main__':
	option = ""; args = ""; bgf_file = ""; ff_file = ""; out_file = ""; solvent_file = ""; size = ""; ff_file = "";
	suffix = ""; nodes = ""; criteria = 0.0; t = 0; margin = "";
	usage = """
Usage: addSolvent.py -b bgfFile -m margin size -o outFile -s solventFile
Make the given structure solvated on water box.

Options:
	-h		print this help message
	-b bgfFile	REQUIRED. The original structure that to be solvated. Only supports BGF format.
	-f ffFile	REQUIRED. CERIUS2 Type Force Field file.
	-s size		OPTIONAL. The size of the solvent box. "x-size y-size z-size"
	-m "n n n"	OPTIONAL. The margin from the molecule to the box in Angstrom units. "x-margin y-margin z-margin"
	-o outFile	REQUIRED. Savename.
	-t solventFile	Optional. BGF box filled with solvent. Default is equilibrated F3C water box by Tod Pascal.
			Removing bad contacts btwn HETATM and the solvent will only be supported for F3C water box.
	"""

	if len(sys.argv) < 2:
		print usage; sys.exit(0)

	options, args = getopt.getopt(sys.argv[1:], 'hb:m:s:o:t:f:', ['help','bgf=','margin=','size=','output=','solvent=','ff='])
	for option, value in options:
		if option in ('-h', '--help'):
			print usage; sys.exit(0)
		elif option in ('-b', '--bgf'):
			bgf_file = value
		elif option in ('-m', '--margin'):
			margin = value
		elif option in ('-s', '--size'):
			size = value
		elif option in ('-o', '--out'):
			out_file = value
		elif option in ('-t', '--solvent'):
			solvent_file = value
		elif option in ('-f', '--ff'):
			ff_file = value
		elif option in (''):
			print usage; sys.exit(0)

	if size != "" and margin != "":
		nu.die("Cannot state size and margin simultaneously. Quitting.")

	# parsing margin
	#margin = (float(margin), float(margin), float(margin))
	if margin != "":
		margin = margin.split(" ")
		margin = [float(margin[0]), float(margin[1]), float(margin[2])]

	if size != "":
		size = size.split(" ")
		size = [float(size[0]), float(size[1]), float(size[2])]

	# setting up defaults
	if out_file == "": out_file = bgf_file[:-4] + "_solv" + bgf_file[-4:]
	if solvent_file == "": solvent_file = "/home/noische/scripts/dat/WAT/f3c_box.bgf"
	if ff_file == "": nu.die("No Force Field File!")

	# main call without silent
	addsolvent(bgf_file, solvent_file, size, margin, out_file, ff_file, False)
