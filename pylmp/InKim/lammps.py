#!/home/noische/program/python27/bin/python
"""
lammps.py
Original: Jan 01 2011 In Kim

Module containing BGF-file information extraction tools including
getAtomBondsID(class BgfFile)
getAtomAnglePairsID(class BgfFile)
getAtomDihedralPairsID(class BgfFile)
"""

import sys
import os
import string
import getopt
from bgf import *
from nutils import *
import dreiding
import pprint
import nutils as nu

version = '110101'

def getAtomFFTypes(myBGF, strip=1):
	"""
	getAtomFFTypes(myBGF, strip=1):
	Returns a list of force field types in the BGF file.
		len(getAtomFFTypes(myBGF)) : the number of force field types in the BGF file.
	"""
	l_fftypes = [];
	for atom in myBGF.a:
		if strip == 0:
			l_fftypes.append(atom.ffType)
		else:
			l_fftypes.append(atom.ffType.strip("- "))

	l_fftypes = removeRepeat(l_fftypes)

	return l_fftypes


def getAtomBondsID(myBGF):
	"""
	getAtomBonds(myBGF):
	Returns a list of pairs of atom numbers (aNo) defined in BgfFile class that is connected each other.
	"""
	# listing atom bonds
	l_pairs = []
	for atom in myBGF.a:
		for element in atom.CONECT:
			a = [atom.aNo, element]
			l_pairs.append(a)

	l_pairs = removeReverse(l_pairs)
	return l_pairs


def getAtomAnglePairsID(myBGF):
	"""
	getAtomAnglePairs(myBGF):
	Returns a list of three pairs of atom numbers (aNo) defined in BgfFile class that is connected each other.
	Three atom pairs are used to add a angle term.
	"""
	# listing atom angles
	l_angles = []
	for atom in myBGF.a:
		for element in atom.CONECT:
			secondatom = myBGF.getAtom(element)
			for element2 in secondatom.CONECT:
				if element2 != atom.aNo:
					a = [atom.aNo, secondatom.aNo, element2]
					l_angles.append(a)

	l_angles = removeReverse(l_angles)
	return l_angles


def getAtomDihedralPairsID(myBGF):
	"""
	getAtomDihedralPairs(myBGF):
	Returns a list of four pairs of atom numbers (aNo) defined in BgfFile class that is connected each other.
	Four atom pairs are used to add a dihedral term.
	"""
	# listing atom dihedrals
	l_diheds = []
	for atom in myBGF.a:
		for element in atom.CONECT: # element contains 2nd atom
			secondatom = myBGF.getAtom(element)
			for element2 in secondatom.CONECT: # from the 2nd atom, (element2 contains the 3rd atom)
				if element2 != atom.aNo: # if 3rd atom is different from 1st atom, (2nd != 3rd)
					thirdatom = myBGF.getAtom(element2) # get the 3rd atom
					for element3 in thirdatom.CONECT:
						fourthatom = myBGF.getAtom(element3)
						if element3 != secondatom.aNo:
							temp_dihed = [atom.aNo, secondatom.aNo, thirdatom.aNo, fourthatom.aNo]
							l_diheds.append(temp_dihed)

	l_diheds = removeReverse(l_diheds)
	return l_diheds


def getAtomBondsFFType(myBGF, strip=1):
	"""
	getAtomBondsFFType(myBGF):
	"""
	l_pairs = getAtomBondsID(myBGF)

	return getFFTypeFromList(myBGF, l_pairs, strip)


def getAtomAnglesFFType(myBGF, strip=1):
	"""
	getAtomAnglesFFType(myBGF):
	"""
	l_angles = getAtomAnglePairsID(myBGF)

	return getFFTypeFromList(myBGF, l_angles, strip)


def getAtomDihedralsFFType(myBGF, strip=1):
	"""
	getAtomDihedralsFFType(myBGF):
	"""
	l_diheds = getAtomDihedralPairsID(myBGF)

	return getFFTypeFromList(myBGF, l_diheds, strip)


def getFFTypeFromList(myBGF, l_pairs, strip=1):
	"""
	getFFTypeFromList(myBGF, l_pairs):
	Returns a list of pairs of force field type.
	- myBGF: a BgfFile class
	- l_pairs: a list of lists that contains a bond, angle, and torsion information in atom number
	- strip: if strip is 1, the string is stripped(trimmed).
	"""
	l_pairs_fftype = [];
	for element in l_pairs:
		l_fftype_small = [];
		for id in element:
			temp = myBGF.getAtom(id)
			if strip == 0:
				str_temp = temp.ffType
			else:
				#str_temp = string.strip(temp.ffType)
				str_temp = temp.ffType.strip("- ")
			l_fftype_small.append(str_temp)
		l_pairs_fftype.append(l_fftype_small)

	l_pairs_fftype = removeReverse(l_pairs_fftype)

	return l_pairs_fftype


def getFFType(myBGF, l_pairs, strip=1):
	"""
	getFFTypeFromList(myBGF, l_pairs):
	Returns a list of pairs of force field type.
	- myBGF: a BgfFile class
	- l_pairs: a list of lists that contains a bond, angle, and torsion information in atom number
	- strip: if strip is 1, the string is stripped(trimmed).
	"""
	l_pairs_fftype = [];
	for element in l_pairs:
		l_fftype_small = [];
		for id in element:
			temp = myBGF.getAtom(id)
			if strip == 0:
				str_temp = temp.ffType
			else:
				#str_temp = string.strip(temp.ffType)
				str_temp = temp.ffType.strip("- ")
			l_fftype_small.append(str_temp)
		l_pairs_fftype.append(l_fftype_small)

	return l_pairs_fftype

def countSimilarLists(list, pattern):
	"""
	countSimilarLists(list, pattern):
	Counts a type of lists from the given list.
	"""
	count = 0;

	for element in list:
		# if element matches to the pattern;
		if compareListPattern(pattern, element):
			count += 1

	return count;


def fixPairCoeff(in_file, out_file, silent=False):
	"""
def LAMMPS_fixPairCoeff():
	This script fixes the bug in the Tod's createLammpsInput.pl script
	that pair_coeff is not fully written when using charmm type pair coefficient
	with dreiding/hbond.

Function Parameters:
	in_file		A filename of Tod's LAMMPS input file
	out_file	A filename of fixed LAMMPS input file
	"""
	# initialize
	d_pair_coeff = dict();
	l_pair_coeff = [];
	d_charmm_pair_coeff = dict();
	l_charmm_pair_coeff = [];
	max_type_num = 0;
	charmm_coeff_type_name = "";

	# open infile
	f = open(in_file)
	while 1:
		line = f.readline()
		if not line:
			break;

		if "pair_coeff" in line:
			### pair_coeff      0:1    1:1    2:lj/charmm/coul/long/opt        3:0.015200000000000         4:2.846421404458384
			parse = line.split()
			parse[1] = int(parse[1])
			parse[2] = int(parse[2])
			l_pair_coeff.append(parse[1:])

			# what is the maximum atom type number?
			if max_type_num < int(parse[1]):
				max_type_num = int(parse[1])
			elif max_type_num < int(parse[2]):
				max_type_num = int(parse[2])
	f.close()

	# get charmm type pair coeffs
	for i in l_pair_coeff:
		if "charmm" in i[2]:
			d_charmm_pair_coeff[i[0], i[1]] = [i[2], float(i[3]), float(i[4])]
	### REMARK: d_charmm_pair_coeff = {(1, 1): ['lj/charmm/coul/long/opt', 0.0152, 2.846421], ...

	# calculate missing terms
	for i in range(1, max_type_num + 1):
		for j in range(i, max_type_num + 1):
			if not d_charmm_pair_coeff.has_key((i, j)):
				epsilon = math.sqrt(d_charmm_pair_coeff[(i, i)][1] * d_charmm_pair_coeff[(j, j)][1]);
				sigma = math.sqrt(d_charmm_pair_coeff[(i, i)][2] * d_charmm_pair_coeff[(j, j)][2]);
				d_charmm_pair_coeff[(i, j)] = [d_charmm_pair_coeff[(i, i)][0], epsilon, sigma]

	# convert dictionary into list
	for i in d_charmm_pair_coeff:
		temp = list(i) + d_charmm_pair_coeff[i]
		temp = nu.flatten(temp)
		l_charmm_pair_coeff.append(temp)

	# update l_pair_coeff with d_charmm_pair_coeff dictionary
	l_new_pair_coeff = [];
	for i in l_pair_coeff:
		if not "charmm" in i[2]:
			l_new_pair_coeff.append(i)
	l_new_pair_coeff = l_new_pair_coeff + l_charmm_pair_coeff
	l_new_pair_coeff.sort()

	##for key, value in sorted(d_pair_coeff.items(), key=lambda item: item[0]):
	##	print key, value

	f1 = open(in_file)
	f2 = open(out_file, 'w')
	flag = False;	# True: pair_coeff is updated, False: not updated

	# save infile
	while 1:
		line = f1.readline()
		if not line:
			break;

		if "pair_coeff" in line:
			if flag:
				pass;
			else:
				# (1, 1, 'lj/charmm/coul/long/opt') [0.0152, 2.846421404458384]
				f2.write("# this pair_coeff is fixed with LAMMPS_fixPairCoeff.py script.\n")
				for i in l_new_pair_coeff:
					line = "pair_coeff\t";
					line += '\t'.join([str(j) for j in i[:3]])
					line += '\t'
					#line += '\t'.join(["{0:>20}".format(j) for j in i[3:]])
					#line += str(key[0]) + "\t" + str(key[1]) + "\t" + str(key[2]) + "\t"
					if "charmm" in i[2]:
						for j in i[3:]:
							line += "{0:>10.6f}".format(j)
					else:
						for k in i[3:]:
							line += "{0:<10}".format(k) + "\t"
					line += "\n"
					f2.write(line)
				flag = True;
		else:
			f2.write(line)
	f1.close()
	f2.close()

	### end of LAMMPS_fixPairCoeff


#-------------------------------------
#
# - Main Function
#
if __name__ == '__main__':

	option = ""; args = ""; bgf_file = ""; mod_file = ""; out_file = ""
	usage = """
Usage: lammps.py -b bgfFile -f forcefield -s suffix
	"""

	options, args = getopt.getopt(sys.argv[1:], 'hb:f:s:', ['help','bgf=','forcefield=','suffix='])
	for option, value in options:
	        if option in ('-h', '--help'):
	                print usage; sys.exit(0)
	        elif option in ('-b', '--bgf'):
	                bgfFile = value
		elif option in ('-f', '--forcefield'):
			ffFile = value
	        elif option in ('-s', '--suffix'):
	                suffix = value
	        elif option in (''):
	                print usage; sys.exit(0)

	##### Initialize
	output = ""
	if suffix == "": suffix = "lammps"


	##### Open a new data file
	myBGF = BgfFile(bgfFile)
	forcefield = dreiding.loadFF(ffFile)
	
	##### Description of the file
	output += "Created by lammps.py\n"
	output += "\n"	


	##### the number of atoms
	output += "{0:8d}".format(len(myBGF.a)) + " atoms" + "\n"


	##### the number of bonds
	output += "{0:8d}".format(len(getAtomBondsID(myBGF))) + " bonds" + "\n"


	##### the number of angles
	output += "{0:8d}".format(len(getAtomAnglePairsID(myBGF))) + " angles" + "\n"


	##### the number of dihedrals
	output += "{0:8d}".format(len(getAtomDihedralPairsID(myBGF))) + " dihedrals" + "\n"


	##### the number of impropers
	### << NOT YET>>
	output += "# Impropers are not yet implemented." + "\n"
	output += "\n"


	##### the number of atom types
	#atomtype = dreiding.loadAtomTypes(forcefield)
	l_fftypes = getAtomFFTypes(myBGF, 1)	# from the BGF file
	output += "{0:8d}".format(len(l_fftypes)) + " atom types" + "\n"


	##### FOR PAIR COEFFS #####
	l_pairfftypes_ff = dreiding.loadPairTypes(forcefield)

	requiredPairFFTypes = [];
	for item in l_pairfftypes_ff:
		if item['LABEL'] in l_fftypes:
			if item not in requiredPairFFTypes:
				requiredPairFFTypes.append(item)
	#print(l_fftypes)
	#print(requiredPairFFTypes)

	##### FOR BOND COEFFS #####
	l_bondfftypes_bgf = getAtomBondsFFType(myBGF, 1)	# NOTE: bonding fftypes read from the BGF file
	l_bondfftypes_ff = dreiding.loadBondTypes(forcefield)	# NOTE: bonding fftypes read from dreiding

	requiredBondFFTypes = [];
	for bgfitem in l_bondfftypes_bgf:
		for ffitem in l_bondfftypes_ff:
			if compareListPattern(ffitem['LABEL'], bgfitem):
				if ffitem not in requiredBondFFTypes:
					requiredBondFFTypes.append(ffitem)
				#if ffitem['ID'] not in bgfitem:
				#	bgfitem.append(ffitem['ID'])
	#print(l_bondfftypes_bgf)
	#print(requiredBondFFTypes)

	##### the number of bond types
	output += "{0:8d}".format(len(requiredBondFFTypes)) + " pair types" + "\n"


	##### the number of angle types
	l_anglefftypes_bgf = getAtomAnglesFFType(myBGF, 1)	# NOTE: fftypes read from the BGF file
	l_anglefftypes_ff = dreiding.loadAngleTypes(forcefield)	# NOTE: fftypes read from dreiding

	requiredAngleFFTypes = [];	# TYPE: list of the dictionary. NOTE: len(requiredAngleFFTypes) := the number of angle types
	for bgfitem in l_anglefftypes_bgf:
		for ffitem in l_anglefftypes_ff:
			if compareListPattern(ffitem['LABEL'], bgfitem):	# if the angle type is found on the forcefield (dreiding),
				if ffitem not in requiredAngleFFTypes:		# reduce the dictionary
					requiredAngleFFTypes.append(ffitem)
				if ffitem['ID'] not in bgfitem:			# reduce the list
					bgfitem.append(ffitem['ID'])		# REMARK: list which the ID for the dictionary is appended
	#print(l_anglefftypes_bgf)
	#print(requiredAngleFFTypes)
	output += "{0:8d}".format(len(requiredAngleFFTypes)) + " angle types" + "\n"
				

	##### the number of dihedral types
	l_dihedralfftypes_bgf = getAtomDihedralsFFType(myBGF, 1)	# NOTE: fftypes read from the BGF file
	l_dihedralfftypes_ff = dreiding.loadTorsionTypes(forcefield)	# NOTE: fftypes read from dreiding

	requiredDihedralFFTypes = [];	# TYPE: list of the dictionary.
	for bgfitem in l_dihedralfftypes_bgf:
		for ffitem in l_dihedralfftypes_ff:
			if compareListPattern(ffitem['LABEL'], bgfitem):
				if ffitem not in requiredDihedralFFTypes:
					requiredDihedralFFTypes.append(ffitem)
				if ffitem['ID'] not in bgfitem:
					bgfitem.append(ffitem['ID'])
	output += "{0:8d}".format(len(requiredDihedralFFTypes)) + " dihedral types" + "\n"
	#print(getAtomDihedralPairsID(myBGF))
	#print(l_dihedralfftypes_bgf)
	#print(requiredDihedralFFTypes)


	##### the number of improper types
	# not yet implemented
	output += "# impropers are not yet implemented.\n"
	output += "\n"


	##### simulation box
	bgfsize = getBGFSize(myBGF, 1.0)	# add 1.0A margin
	output += "\t" + "{0:10.6f}".format(bgfsize[0]) + "  " + "{0:10.6f}".format(bgfsize[1]) + " xlo xhi" + "\n"
	output += "\t" + "{0:10.6f}".format(bgfsize[2]) + "  " + "{0:10.6f}".format(bgfsize[3]) + " ylo yhi" + "\n"
	output += "\t" + "{0:10.6f}".format(bgfsize[4]) + "  " + "{0:10.6f}".format(bgfsize[5]) + " zlo zhi" + "\n"
	output += "\n"


	##### Masses
	#l_fftypes = getAtomFFTypes(myBGF, 1)
	atomtypes = dreiding.loadAtomTypes(forcefield)

	output += "Masses" + "\n"
	for index, item in enumerate(l_fftypes):
		output += "{0:8d}".format(index + 1) + "\t" + "{0:>10.5f}".format(atomtypes[item]['MASS']) + "\t\t" + "# " + str(item) + "\n"	# REMARK: id = index + 1
	output += "\n"


	##### Pair Coeffs
	output += "Pair Coeffs" + "\n"
	for index, item in enumerate(l_fftypes):
		for ffitem in requiredPairFFTypes:
			if item == ffitem['LABEL']:
				ffitem['ID'] = index + 1
				output += "{0:8d}".format(ffitem['ID']) + "\t" + "{0:10.6f}".format(ffitem['EPSILON']) + "\t" + "{0:10.6f}".format(ffitem['SIGMA']) + "\t\t" + "# " + ffitem['LABEL'] + "\n"
	output += "\n"

	#print(requiredPairFFTypes)

	##### Bond Coeffs
	output += "Bond Coeffs" + "\n"
	for index, item in enumerate(l_bondfftypes_bgf):
		for ffitem in requiredBondFFTypes:
			#if item == ffitem['LABEL'] or rev_item == ffitem['LABEL']:
			if compareListPattern(ffitem['LABEL'], item):
				ffitem['ID'] = index + 1
				output += "{0:8d}".format(ffitem['ID']) + "\t" + "{0:10.6f}".format(ffitem['K']) + "\t" + "{0:10.6f}".format(ffitem['R']) + "\t\t" + "# " + str(ffitem['LABEL']) + "\n"
	output += "\n"

	##### Angle Coeffs
	output += "Angle Coeffs" + "\n"
	for index, item in enumerate(requiredAngleFFTypes):
		item['ID'] = index + 1
		output += "{0:8d}".format(item['ID']) + "\t" + "{0:10.6f}".format(item['K']) + "\t" + "{0:10.6f}".format(item['THETA']) + "\t\t" + "# " + str(item['LABEL']) + "\n"	# REMARK: id = index + 1
	output += "\n"

	##### Dihedral Coeffs
	output += "Dihedral Coeffs" + "\n"

	wholeDihedralFFTypes = [];
	for item in getAtomDihedralPairsID(myBGF):
		type = getFFTypeFromList(myBGF, [item])
		wholeDihedralFFTypes.append(type[0])

	for index, item in enumerate(requiredDihedralFFTypes):
		item['ID'] = index + 1
		item['OCCURRENCE'] = countSimilarLists(wholeDihedralFFTypes, item['LABEL'])
		output += "{0:8d}".format(item['ID']) + "\t" + "{0:10.6f}".format(item['K']/float(item['OCCURRENCE'])/2.0) + "\t" + "{0:10d}".format(int(item['n'])) + "\t" + "{0:10d}".format(int(0)) + "\t\t" + "# " + str(item['LABEL']) + "\n"	# charmm dihedral style

	output += "\n"
	
	#print(requiredDihedralFFTypes)

	##### Atoms
	## output: atom-ID molecule-ID atom-type q x y z
	output += "Atoms" + "\n"
	for item in myBGF.a:
		output += "{0:8d}".format(item.aNo) + "{0:8d}".format(item.rNo) + "{0:8d}".format(l_fftypes.index(str.strip(item.ffType))+1) 
		output += "{0:12.5f}".format(item.charge) + "{0:12.5f}".format(item.x) + "{0:12.5f}".format(item.y) + "{0:12.5f}".format(item.z)
		output += "\t" + "# " + item.ffType + "\n"
	output += "\n"

	##### Bonds
	## output: ID type atom1 atom2
	output += "Bonds" + "\n"
	#output += "#	ID	type	atom1	atom2" + "\n"
	for index, item in enumerate(getAtomBondsID(myBGF)):
		id = index + 1
		type = getFFTypeFromList(myBGF, [item])		## REMARK: requires improvement!!
		type_id = 0;
		for ffitem in requiredBondFFTypes:
			if compareListPattern(ffitem['LABEL'], type[0]):
				type_id = ffitem['ID']
		output += "{0:8d}".format(id) + "{0:8d}".format(type_id) + "{0:8d}".format(item[0]) + "{0:8d}".format(item[1]) + "\t" + "# " + str(type[0]) + "\n"
	output += "\n"

	##### Angles
	## output: ID type atom1 atom2 atom3
	output += "Angles" + "\n"
	for index, item in enumerate(getAtomAnglePairsID(myBGF)):
		id = index + 1
		type = getFFTypeFromList(myBGF, [item])		## REMARK: requires improvement!!
		type_id = 0;
		for ffitem in requiredAngleFFTypes:
			if compareListPattern(ffitem['LABEL'], type[0]):
				type_id = ffitem['ID']
		output += "{0:8d}".format(id) + "{0:8d}".format(type_id) + "{0:8d}".format(item[0]) + "{0:8d}".format(item[1]) + "{0:8d}".format(item[2]) + "\t" + "# " + str(type[0]) + "\n"
        output += "\n"

	##### Dihedrals
	## output: ID type atom1 atom2 atom3 atom4
	output += "Dihedrals" + "\n"

	for index, item in enumerate(getAtomDihedralPairsID(myBGF)):
		id = index + 1
		type = getFFTypeFromList(myBGF, [item])		## REMARK: requires improvement!!
		type_id = 0;
		for ffitem in requiredDihedralFFTypes:
			if compareListPattern(ffitem['LABEL'], type[0]):
				type_id = ffitem['ID']

		output += "{0:8d}".format(id) + "{0:8d}".format(type_id) + "{0:8d}".format(item[0]) + "{0:8d}".format(item[1]) + "{0:8d}".format(item[2]) + "{0:8d}".format(item[3]) + "\t" + "# " + str(type[0]) + "\n"
        output += "\n"

	##### Impropers
	### not yet implemented
	#print(output)

	file = open("data." + suffix, "w")
	file.write(output)
	file.close()

	##### control file
	output = ""
	output += "units" + "\t\t" + "real" + "\n"
	output += "atom_style" + "\t" + "full" + "\n"

	##### determine the boundary condition
	bc = ""
	for letter in myBGF.PERIOD:
		if letter == "0": bc += "f "
		elif letter =="1": bc += "p "
	if myBGF.PERIOD == "": bc = "f f f"
	output += "boundary" + "\t" + bc + "\n"

	output += "dielectric" + "\t" + "1" + "\n"
	output += "special_bonds" + "\t" + "lj/coul 0.0 0.0 1.0" + "\n"
	output += "\n"

	output += "pair_style" + "\t"
	if not myBGF.PERIOD:	# finite
		output += "lj/charmm/coul/charmm 9.0 10.0" + "\n"
	else:
		output += "lj/charmm/coul/long/opt 9.0 10.0" + "\n"
	output += "bond_style" + "\t" + "harmonic" + "\n"
	output += "angle_style" + "\t" + "harmonic" + "\n"
	output += "dihedral_style" + "\t" + "charmm" + "\n"
	output += "improper_style" + "\t" + "none" + "\n"	## REMARK: requires revision
	output += "kspace_style" + "\t"

	if not myBGF.PERIOD:	# finite
		output += "none" + "\n"
	else:
		output += "pppm 0.0001" + "\n"

	output += "\n"
	output += "read_data" + "\t" + "data." + suffix + "\n"

	output += "\n"
	output += "pair_modify" + "\t" + "mix geometric" + "\n"
	output += "neighbor" + "\t" + "2.0 multi" + "\n"
	output += "neigh_modify" + "every 2 delay 4 check yes" + "\n"
	output += "thermo_style" + "\t" + "multi"

	output += "\n"
	output += "variable" + "\t" + "input index in." + suffix + "\n"
	output += "variable" + "\t" + "sname index " + suffix + "\n"

	file = open("in." + suffix, "w")
	file.write(output)
	file.close()
	#print(output)
