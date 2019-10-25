#!/home/noische/program/python27/bin/python

import sys
import re
import string
import getopt
import optparse
import math
import time
from os import popen

# Pizza.py module
sys.path.append("/home/noische/program/pizza-1Oct10")
import dump

option = ""; args = ""; bgf_file = ""; trj_file = ""; out_file = ""
usage = """
getLAMMPSRg.py: 
	Read atom data from the LAMMPS trajectory file and atom mass from LAMMPS data file
	Write the Radius of gyration to the output file

Usage: getLAMMPSRg.py -d dat_file -t trj_file -a atom_selection -o out_file

	dat_file: REQUIRED. LAMMPS data file for atom mass information.
	trj_file: REQUIRED. LAMMPS trajectory file for Rg calculation
	atom_selection: REQUIRED. "a-b": atom numbers within a and b. Other selections are not supported.
	out_file: OPTIONAL. Output file. (Default: suffix_Rg.dat)

	e.g.) getLAMMPSRg.py -d data.p16_02 -a "1-2344" -t p16_02.heat.lammpstrj -o p16_02_Rg.dat
		will calculate Rg of atom id 1 to 2344.
"""

version = "111031"

def getLAMMPSRg(trj_file, dat_file, out_file, atoms, silent=True):
	"""
getLAMMPSDensity: read atom data from the LAMMPS trajectory file and atom mass from LAMMPS data file
           write the density to the output file
	"""
	if not silent: print("Getting the density from the trajectory " + trj_file + " with the information in " + dat_file + " to " + out_file)

	# initialization
	mass = 0; l_information = []; 

	# get masses from data file
	atom_mass = dict();
	atom_type_number = 0;
	f_dat_file = open(dat_file)
	while 1:
		line = f_dat_file.readline()

		if not line:
			break;

		# determine how many atoms in the data file
		if "atom type" in line:
			parse = re.split('\s*', line)
			atom_type_number = int(parse[1])

		# get masses of atom type as a dictionary
		if "Masses" in line:
			line2 = f_dat_file.readline()
			for i in range(0, atom_type_number):
				line2 = f_dat_file.readline()
				parse = re.split('\s*',line2)

				atom_mass[int(parse[1])] = float(parse[2])
	if not silent: print("Atom masses loaded from the LAMMPS data file " + dat_file)

	# parse the atom selection
	selected_atoms = [];
	[a, b] = [int(i) for i in atoms.split("-")]
	selected_atoms = range(a, b+1)
	if not silent: print(str(len(selected_atoms)) + " atoms are selected. (from atom id " + str(a) + " ~ " + str(b) + ")")

	# get volumes and atoms at the specific step from the trajectory

	myDUMP = dump.dump(trj_file, 0)
	while 1:
		time = myDUMP.next()
		if time == -1:
			sys.stdout.flush()
			if not silent: print("\nDone!")
			break;

		sys.stdout.write("\rReading timestep.. " + str(time))
		sys.stdout.flush()

		l_timestep = myDUMP.time()
		atoms = myDUMP.viz(len(l_timestep) - 1)[2]

		# delete other atom info
		n_atoms = len(atoms)
		rg_atoms = []
		for index, atom in enumerate(atoms):
			if int(atom[0]) >= a and int(atom[0]) <= b:
				rg_atoms.append(atom)
		del(atoms)
		atoms = rg_atoms
		del(rg_atoms)

		# mass
		mass = 0;
		for atom in atoms:
			mass += atom_mass[atom[1]]

		# center of mass
		m_x = 0; m_y = 0; m_z = 0;
		for atom in atoms:
			m_x += atom[2] * atom_mass[atom[1]] / mass
			m_y += atom[3] * atom_mass[atom[1]] / mass
			m_z += atom[4] * atom_mass[atom[1]] / mass

		# radius of gyration
		rg = 0;
		for atom in atoms:
			r = (atom[2] - m_x)**2 + (atom[3] - m_y)**2 + (atom[4] - m_z)**2
			rg += r * atom_mass[atom[1]] / mass
		rg = math.sqrt(rg)
		l_information.append([time, rg])

	# file output
	f_out_file = open(out_file, "w")
	f_out_file.write("Timestep\tRg\n")
	for i in l_information:
		f_out_file.write("{0:>10}\t{1:>10.4f}\n".format(*i))
	f_out_file.write("="*60 + "\n")
	f_out_file.close()

	if not silent: print("SUCCESS: The Rg information is written in " + out_file)

	return 1;	# success

	### end of updatebgf


if __name__ == "__main__":

	if len(sys.argv) < 2:
		print(usage);
		sys.exit(1)

	dat_file = ""; trj_file = ""; atoms = ""; out_file = "";

	options, args = getopt.getopt(sys.argv[1:], 'hd:t:o:a:', ['help','data=','trj=','out=','atom='])
	for option, value in options:
		if option in ('-h', '--help'):
			print(usage)
			sys.exit(0);
		elif option in ('-d', '--data'):
			dat_file = value
		elif option in ('-t', '--trj'):
			trj_file = value
		elif option in ('-o', '--out'):
			out_file = value
		elif option in ('-a', '--atom'):
			atoms = value
		elif option == NULL:
			print(usage)
			sys.exit(0)

	if trj_file == "":
		print("You have to specify trajectory file.")
		sys.exit(1)
	if dat_file == "":
		print("You have to specify data file for atom mass.")
		sys.exit(1)
	if out_file == "":
		out_file = trj_file[5:] + "_Rg.dat"

	# main call
	print("\n" + sys.argv[0] + " version " + str(version) + "\n")
	getLAMMPSRg(trj_file, dat_file, out_file, atoms, silent=False)


