#!/home/noische/program/python27/bin/python
"""
patch_crosslink.py
Original: Jun 02 2011 In Kim

"""

# Python Modules
import sys
import os
import string
import shutil
import time
import copy
import lammps

# Globals
version = "110602"

def pair_coeff_patch_wo_xlinker(in_file, silent):
	"""
pair_coeff_patch_wo_xlinker():
	Patches the LAMMPS input file created by Tod's script on qch.
	This is because there is a problem on Tod's createLammpsInput.pl script which sometimes omits
	describing lj/charmm/coul/long/opt with hb/dreiding/lj pair_coeff term.
	"""

	out_file = in_file + "_patch"

	# pair coeff patch
	lammps.fixPairCoeff(in_file, out_file)
	shutil.copy(out_file, in_file)

	f_in_file = open(in_file)
	f_out_file = open(out_file, 'w')

	while 1:
		line = f_in_file.readline()
		if not line:
			break;

		elif "hbond/dreiding/lj 2 5.0 90" in line:
			f_out_file.write("pair_style      hybrid/overlay hbond/dreiding/lj 2 4.5 5.0 90 lj/charmm/coul/long/opt  9.0 10.0\n")
		elif "#$#" in line:
			line = line.replace("#$#", "")
			f_out_file.write(line)
		#elif "kspace_style    pppm 0.0001" in line:
		#	f_out_file.write("kspace_style    pppm 0.001\n")
		elif "image yes" in line:
			pass;
		elif "1 all atom 25" in line:
			f_out_file.write("dump            1 all custom 25 ${sname}_min.lammpstrj id type xu yu zu ix iy iz\n")
		else:
			f_out_file.write(line)

	f_in_file.close()
	f_out_file.close()

	# data file patch
	datafile = "data" + in_file[2:]
	os.system("sed -i 's/0 # X/0 0 # X/' " + datafile)
	os.system("sed -i 's/Impropers//' " + datafile)

	shutil.copy(out_file, in_file)
	if not silent: print("pair_coeff_patch (version " + version + ") is applied on " + in_file)


def pair_coeff_patch_w_xlinker(in_file, silent):
	"""
pair_coeff_patch_w_xlinker():
	Patches the LAMMPS input file created by Tod's script on qch.
	This is because there is a problem on Tod's createLammpsInput.pl script which sometimes omits
	describing lj/charmm/coul/long/opt with hb/dreiding/lj pair_coeff term.
	"""

	out_file = in_file + "_patch"

	# pair coeff patch
	lammps.fixPairCoeff(in_file, out_file)
	shutil.copy(out_file, in_file)

	f_in_file = open(in_file)
	f_out_file = open(out_file, 'w')

	while 1:
		line = f_in_file.readline()
		if not line:
			break;

		elif "hbond/dreiding/lj 2 5.0 90" in line:
			f_out_file.write("pair_style      hybrid/overlay hbond/dreiding/lj 2 4.5 5.0 90 lj/charmm/coul/long/opt  9.0 10.0\n")
		elif "kspace_style    pppm 0.0001" in line:
			f_out_file.write("kspace_style    pppm 0.001\n")
		elif "image yes" in line:
			pass;
		elif "1 all atom 25" in line:
			f_out_file.write("dump            1 all custom 25 ${sname}_min.lammpstrj id type xu yu zu ix iy iz\n")
		else:
			f_out_file.write(line)

	f_in_file.close()
	f_out_file.close()

	# data file patch
	datafile = "data" + in_file[2:]
	os.system("sed -i 's/0 # X/0 0 # X/' " + datafile)
	os.system("sed -i 's/Impropers//' " + datafile)

	shutil.copy(out_file, in_file)
	if not silent: print("pair_coeff_patch (version " + version + ") is applied on " + in_file)
	

if __name__ == "__main__":

	usage = """
patch_crosslink.py: Do not access this file from the command line.
	"""

	print(usage)
	sys.exit(0)

