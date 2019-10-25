#!/home/noische/program/python27/bin/python
"""
template.py
Original: Dec 28 2011 In Kim
"""

# Python Modules
import sys
import os
import string
import random
import time
import getopt
import pprint

# Custom Modules
sys.path.append("/home/noische/scripts")
sys.path.append("/home/noische/script")
import bgf
import bgftools
import nutils as nu
import scipy

# Globals
version = '111228'

def getViscosity(log_file, profile_file, out_file, silent=False):
	"""
def template():
	Write what this function does.

Function Parameters:
	log_file	A string of filename or BgfFile class.
	profile_file	A string of filename or BgfFile class.
	"""
	# Initialize
	log_data = [];
	profile_data = [];
	boxlength = 0;
	
	f_out_file = open(out_file, 'w')
	

	# Load log_file (dP)
	f_log_file = open(log_file);
	while 1:
		line = f_log_file.readline()
		if not line:
			break;
		if "Step dp TotEng Temp" in line:
			break;
		if "Box length" in line and "print" not in line:
			parse = line.split()
			boxlength = float(parse[-1])

	while 1:
		line = f_log_file.readline()
		if not line:
			break;

		# log_data: [Step dp TotEng Temp]
		parse = line.split()
		if len(parse) != 4:
			break;

		log_data.append([int(parse[0]), float(parse[1]), float(parse[2]), float(parse[3])])


	# Load .profile and calculate dvx/dvz
	f_profile_file = open(profile_file);

	while 1:
		timestep = 0; bin = 0;
		vx = []; vz = [];

		line = f_profile_file.readline()
		if not line:
			break;

		if "#" in line:
			continue;

		parse = line.split()
		if len(parse) == 2:
			timestep = int(parse[0])
			bin = int(parse[1])
			
		# read vz-vx pairs
		for i in range(0, bin):
			dataline = f_profile_file.readline()
			parsedata = dataline.split()
			vz.append(float(parsedata[1])*boxlength)
			vx.append(float(parsedata[3]))

		if len(vz) != bin or len(vx) != bin:
			nu.die("The number of vectors for linear regression in the profile file does not match.")

		# regression of vx wrt vz (i.e. x = vz, y = vx in common case)
		(ar, br) = scipy.polyfit(vz, vx, 1)

		temp = [timestep, ar]
		profile_data.append(temp)
		#f_out_file.write(str(temp)+"\n")	# profile reader looks good 2012.2.2

	# merge two data: log and profile
	# merged_data: [Step dp TotEng Temp (dvx/dvz)]
	merged_data = [];
	for item1 in log_data:
		for item2 in profile_data:
			if item1[0] == item2[0]:
				temp = item1 + item2[1:]
				merged_data.append(temp)

	# viscosity = - dp / (dvx/dvz)
	for i in merged_data:
		vis = -1 * i[1] / i[4]
		i.append(vis)

	for i in merged_data:
		line = "";
		for j in i:
			line += str(j) + "\t"
		line += "\n"
		f_out_file.write(line)

	# close files
	f_out_file.close()

	### end of template


if __name__ == '__main__':

	option = ""; args = ""; log_file = ""; size = 0.0; profile_file = ""; out_file = "";
	number = 0
	usage = """
Usage: LAMMPS_getViscosity.py -l logfile -p profile -o output

Options are:
	-b	Input BGF file.
	-o	Output BGF file.
	"""

	if len(sys.argv) < 2:
		print(usage); sys.exit(0)

	options, args = getopt.getopt(sys.argv[1:], 'hl:p:o:', ['help','log=','profile=','out='])
	for option, value in options:
		if option in ('-h', '--help'):
			print usage; sys.exit(0)
		elif option in ('-l', '--log'):
			log_file = value
		elif option in ('-o', '--output'):
			out_file= value
		elif option in ('-p', '--profile'):
			profile_file = value
		elif option in (''):
			print(usage); sys.exit(0)

	# default settings
	if not out_file: out_file = os.path.basename(log_file).split(".log")[0] + "" + ".output"

	getViscosity(log_file, profile_file, out_file, silent=False)

