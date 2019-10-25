#!/home/noische/program/python27/bin/python
"""
PEI_anneal.py
Original: Feb 10 2012

Abstract:

For all BGF file in the directory

Perform 
"""

# Python Modules
import sys
import os
import string
import shutil
import random
import time
import copy
import getopt
import optparse
from types import *

# BGF Modules
sys.path.append("/home/noische/script")
sys.path.append("/home/noische/scripts")
from bgf import *
import bgftools
import bgfselect
import nutils as nu
import LAMMPS_getLogData
import LAMMPS_trj2bgf

# Pizza.py Modules
sys.path.append("/home/noische/program/pizza-1Oct10")
from dump import *

# Cluster settings
from clusterSetting import *

# Globals
version = '111202'

def anneal(bgf_dir, ff_file, suffix, out_file):

	# *** variables ***

	# used for log file output
	output = ""

       	# initialization
	#t1 = 0; t2 = 0;	# clock for measuring time for a cycle
	
	curr_dir = os.path.abspath(bgf_dir)
	temp_dir = curr_dir + "/anneal/"

	# LAMMPS executive command: defined in ~/scripts/clusterSettings.py
	#lammps_parallel = mpi_command + " -np 8 " + lammps_command
	lammps_parallel = "/home/tpascal/codes/openmpi/1.4.3/bin/mpirun -np 8 " + lammps_command

	# Annealing procedure
	anneal_proc = "/home/noische/script/dat/in.lammps.anneal"

	# get BGF file list in the directory
	pei_files = glob.glob(bgf_dir + "/*.bgf")
	pei_files.sort()

	os.chdir(bgf_dir)

	#print("pei_files: ", pei_files)

	# for file in the directory
	for file in pei_files:
		### working filename
		print("*** Working on %s" % file)

		filename = file.split(".bgf")[0]
		print("filename: " + filename)

		# centerBGF
		print("Centering %s" % file)
		cmd_centerbgf = "/home/tpascal/scripts/centerBGF.pl -b " + file + " -f " + ff_file + " -s " + filename + ".center.bgf" + " -c com_center"
		nu.shutup(); os.system(cmd_centerbgf); nu.say();

		filename = filename.split("/")[-1]
		print("filename: " + filename)

		# createLammps: with finite options
		# Note that the script refers /home/noische/script/data/in.lammps.anneal for the annealing process
		cmd_createLammpsInput = "/home/tpascal/scripts/createLammpsInput.pl -b " + filename + ".center.bgf" + " -f " + ff_file + " -t /home/noische/script/data/in.lammps.anneal " + " -s " + filename + " -o finite "
		nu.shutup(); os.system(cmd_createLammpsInput); nu.say();

		# patch in.lammps file
		print("Patching infile")
		cmd_sed = "sed -i 's/2 5.0 90/2 4.5 5.0 90/' in." + filename
		nu.shutup(); os.system(cmd_sed); nu.say();

		# patch data.lammps file
		print("Patching datafile")
		cmd_mv = "mv data." + filename + " data." + filename + ".bak"
		os.system(cmd_mv)
		file1 = open("data." + filename, "w")
		file2 = open("data." + filename + ".bak")
		while 1:
			line = file2.readline()
			if not line: break;
			if "xlo" in line:
				line = "-50.0 50.0 xlo xhi\n"
				file1.write(line)
			elif "ylo" in line:
				line = "-50.0 50.0 ylo yhi\n"
				file1.write(line)
			elif "zlo" in line:
				line = "-50.0 50.0 zlo zhi\n"
				file1.write(line)
			elif "# X N_3 C_3 X" in line:
				line = line.replace("# X N_3 C_3 X", "0 # X N_3 C_3 X")
				file1.write(line)
			elif "# X C_3 C_3 X" in line:
				line = line.replace("# X C_3 C_3 X", "0 # X C_3 C_3 X")
				file1.write(line)
			else:
				file1.write(line)

		# run
		print("Running simulation")
		lammps_thermo_log_file = filename + ".anneal.log"
		cmd_lammps = lammps_parallel + " -log " + lammps_thermo_log_file
		nu.shutup(); os.system(cmd_lammps); nu.say();

		# read potential energy from file
		# in in.lammps.anneal: "poteng" is the averaged value of potential energy after cooldown
		print("Treating Potential Energy")
		poteng = [];
		poteng_files = glob.glob(filename + "*poteng")
		for poteng_file in poteng_files:
			f_poteng_file = open(poteng_file)
			while 1:
				line = f_poteng_file.readline()
				if not line:
					break;
				if not line[0] == "#":
					parse = line.split()
					poteng.append([parse[0], float(parse[1])])
		poteng = nu.removeRepeat(poteng)

		### print potential energy
		for i in poteng:
			print(i)

		# when is the minimum energy during annealing?
		min_poteng = 0; min_poteng_step = 0;	# usually poteng < 0
		for i in poteng:
			if i[1] < min_poteng:
				min_poteng_step = i[0]
				min_poteng = i[1]

		# convert that timestep into BGF
		trj_file = filename + ".lammpstrj"
		anneal_output_bgf_file = filename + ".annealed.bgf"
		result = LAMMPS_trj2bgf.getLAMMPSTrajectory(pei_file, trj_file, anneal_output_bgf_file, min_poteng_step, False)
		print("The lowest potential energy structure found")

		# create LAMMPS input for annealed BGF file

		# run the NVT simulation at 300 K

		# average the potential energy

		# average Rg

		# write data: filename, the minimum structure, avr_poteng, avr_Rg


	"""
	l_init_amines = bgftools.getAmineGroupInfo(initial_bgf_file)
	setAmineGroupInfo = bgftools.setAmineGroupInfo(initial_bgf_file, initial_bgf_file)

	# scratch setting and file copy
	if not os.path.isdir(temp_dir):
		os.makedirs(temp_dir)

	# managing hostfile: "COMMENTED"
	# Important Note: assume that "cat $PBS_NODEFILE > hostfile" is already executed in the .pbs file
	#hostfile = open(temp_dir + "/hostfile", 'w')
	#for i in range(0, 12): hostfile.write(os.environ["HOST"] + "\n")
	#hostfile.close()
	shutil.copy(initial_bgf_file, temp_dir)
	shutil.copy(ff_file, temp_dir)
	os.chdir(temp_dir)

	# LAMMPS input script files for crosslinking
	#init_script_original_file = "/home/noische/research/dendrimer/simulation/600/procedures/in.noinit"
	#init_script_original_file = "/home/noische/research/dendrimer/simulation/600/procedures/in.initialize"
	#crosslink_script_original_file = "/home/noische/research/dendrimer/simulation/600/procedures/in.xlnk_npt_1ps_k16_r28"
	#crosslink_script_original_file = "/home/noische/research/dendrimer/simulation/600/procedures/in.xlnk_npt_1ps_k8_r28"
	#crosslink_script_original_file = "/home/noische/research/dendrimer/simulation/600/procedures/in.xlnk_npt_1ps_k16_r35"

	shutil.copy(init_script_original_file, temp_dir)
	shutil.copy(crosslink_script_original_file, temp_dir)

	init_script_filename = init_script_original_file.split("/")[-1]	# initial script filename
	crosslink_script_filename = crosslink_script_original_file.split("/")[-1]	# crosslinking script filename

	createLammpsInput = "createLammpsInput.pl -b " + bgf_file + " -f " + ff_file + " -s " + suffix + " -o 'no shake' -t " + init_script_filename + " > /dev/null"
	os.system(createLammpsInput)
	initial_in_file = "in." + suffix
	"""



if __name__ == '__main__':

	option = ""; args = ""; bgf_dir = ""; ff_file = ""; out_file = "";
	suffix = ""; nodes = ""; t = 0; ratio = 0.0;
	usage = """
Usage: PEI_anneal.py -d bgfDirectory -f forcefield -s suffix -o log_file
	"""

	if len(sys.argv) < 2:
		print(usage);
		sys.exit(0)

	options, args = getopt.getopt(sys.argv[1:], 'hd:f:s:o:', ['help','directory=','forcefield=','suffix=','output='])
	for option, value in options:
		if option in ('-h', '--help'):
			print usage; sys.exit(0)
		elif option in ('-d', '--directory'):
			bgf_dir = value
		elif option in ('-f', '--forcefield'):
			ff_file = value
		elif option in ('-s', '--suffix'):
			suffix = value
		elif option in ('-o', '--option'):
			out_file = value
		elif option in (''):
			print usage; sys.exit(0)

		if suffix == "": suffix = "anneal"
		if out_file == "": out_file = suffix + "_log.out"
		if ff_file == "": ff_file = "/home/noische/ff/dreiding-den.par"

	anneal(bgf_dir, ff_file, suffix, out_file)


"""
qch.kaist.ac.kr:/disk2/scratch/noische/_simulation/600/crosslink.2/xlnk02001/600xlnk02001f3c.force
can see fix recenter here
"""
