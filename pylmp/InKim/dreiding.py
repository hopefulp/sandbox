#!/home/noische/program/python27/bin/python
"""
dreiding.py
Original: Jan 01 2011 In Kim

Module containing general data of atoms
"""

import os, sys, re, string, pprint
import nutils as nu

version = '110101'

def loadFF(file):

	if file == "":
		#nu.warn("Forcefield type is not specified. Using DREIDING 2 as default.")
		file = "/home/noische/ff/DREIDING2.21.ff"

	file = file.strip()

	lines = []
	try:
		ffFile = open(file)
		while 1:
			line = ffFile.readline()
			if line == "": break
			line = line.strip()
			line = line.rstrip('\n')
			line = re.split('\s*', line)
			lines.append(line)
	except IOError:
		nu.die("Force Field File " + file + " open failed")
	else:
		pass;
		#print("Force Field Dreiding loaded")

	start_index = lines.index(["VERSION"])
	stop_index = lines.index(["END"], start_index)
	temp = nu.flatten(lines[start_index + 1 : stop_index])

	for i in temp:
		if "CERIUS" in i:
			return lines

	nu.die("Forcefield file " + file + " is not a CERIUS2 type file.")


def loadAtomTypes(lines):

	temp_list = [];
	atom_types = dict();

	start_index = lines.index(["ATOMTYPES"])
	stop_index = lines.index(["END"], start_index)
	temp_list = lines[start_index + 1 : stop_index]

	for item in temp_list:
		element = dict();
		label = item[0]
		element['LABEL'] = label
		element['ATOM'] = item[1]
		element['MASS'] = float(item[2])
		element['CHARGE'] = float(item[3])
		element['NUMBONDS'] = int(item[4])	# something is weird.. why numbonds of C is 3, not 4?
		element['OTHER'] = item[5]	# what is this?
		element['LONEPAIRS'] = int(item[6])	# what is this, too?

		atom_types[label] = element	# atom types is the key for the dictionary 'atomtypes'

	return atom_types


def loadPairTypes(lines):

	temp_list = [];
	pair_types = [];

	start_index = lines.index(["DIAGONAL_VDW"])
	stop_index = lines.index(["END"], start_index)
	temp_list = lines[start_index + 1 : stop_index]

	for index, item in enumerate(temp_list):
		element = dict();
		label = string.strip(item[0])
		element['LABEL'] = label
		element['VDWTYPE'] = item[1]
		if str(item[1]) == "LJ_6_12":
			element['SIGMA'] = float(float(item[2]) / 2**(1.0/6.0))	# defined in LAMMPS: pair_style: lj/charmm/coul/long/opt
		else:
			element['SIGMA'] = float(item[2])	# distance units
		element['EPSILON'] = float(item[3])	# energy units
		element['ID'] = index

		pair_types.append(element)

	return pair_types


def loadBondTypes(lines):

	temp_list = [];
	bond_types = [];

	start_index = lines.index(["BOND_STRETCH"])
	stop_index = lines.index(["END"], start_index)
	temp_list = lines[start_index + 1 : stop_index]

	for index, item in enumerate(temp_list):
		element = dict();
		label = [ string.strip(item[0]), string.strip(item[1]) ]
		element['LABEL'] = label
		element['STYLE'] = item[2]
		if item[2] == "HARMONIC":
			element['K'] = float(item[3])/2	# K & R0 is defined in LAMMPS: bond_style harmonic.
			element['R'] = float(item[4])
		element['ID'] = index
		
		bond_types.append(element)

	return bond_types


def loadAngleTypes(lines):

	temp_list = [];
	angle_types = [];

	start_index = lines.index(["ANGLE_BEND"])
	stop_index = lines.index(["END"], start_index)
	temp_list = lines[start_index + 1 : stop_index]

	for index, item in enumerate(temp_list):
		element = dict();
		label = [ string.strip(item[0]), string.strip(item[1]), string.strip(item[2]) ]
		element['LABEL'] = label
		element['STYLE'] = string.strip(item[3])
		element['K'] = float(item[4])/2
		element['THETA'] = float(item[5])
		element['ID'] = index

		angle_types.append(element)

	return angle_types	# type: list


def loadTorsionTypes(lines):

	temp_list = [];
	torsion_types = [];

	start_index = lines.index(["TORSIONS"])
	stop_index = lines.index(["END"], start_index)
	temp_list = lines[start_index + 1 : stop_index]

	for index, item in enumerate(temp_list):
		element = dict();
		label = [ string.strip(item[0]), string.strip(item[1]), string.strip(item[2]), string.strip(item[3]) ]
		element['LABEL'] = label
		element['STYLE'] = string.strip(item[4])
		element['K'] = float(item[5])
		element['n'] = float(item[6])
		element['d'] = float(item[7])
		element['ID'] = index

		torsion_types.append(element)

	return torsion_types	# TYPE: list

##### impropers????

##### Usage #####

#a = loadFF("")
#atomtypes = loadAtomTypes(a)
#pairtypes = loadPairTypes(a)
#bondtypes = loadBondTypes(a)
#angletypes = loadAngleTypes(a)
#torsiontypes = loadTorsionTypes(a)
#pprint.pprint(torsiontypes)
#print(at['N_3'])
