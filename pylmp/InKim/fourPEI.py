#!/home/noische/Enthought/Canopy_64bit/User/bin/python
"""
fourPEI.py
Original: Mar 26 2011 In Kim

# Updates
- 111003:
	* crosslink() function will NOT be called in this code any more: related codes are deleted or commented.
	  fourPEI.py will only perform merging, not performing crosslinking.
- 110514:
	Updated 27 cubic-style crosslinking
- 110527:
	Updated yes_crosslink, silent, Rg, and logging (log_file)
- 110602:
	Merged monomers will be recorded on REMARK.
	Residue numbers of the monomers will lose their original number and have new numbers 1 to 4.
"""

# Python Modules
import sys
import os
import string
import shutil
import random
import time
import glob
import getopt

# Custom Modules
sys.path.append("/home/noische/script")
sys.path.append("/home/noische/scripts")
import bgf
import nutils as nu
import centerBGF
#import crosslink


# Globals
version = '110527'

def fourPEI(bgf_file, ff_file, suffix, out_file, monomers, Rg, directory, silent):
	"""
def fourPEI(bgf_file, ff_file, suffix, out_file, monomers, Rg, directory, silent):
	Runs a simulation which performs a crosslinking with given monomers.

Function Parameters:
	bgf_file	A string which contains a PEI information in a BGF format. 
			(e.g. "file1.bgf file2.bgf file3.bgf ...")

	ff_file		A string which contains a ForceField information. 
			(e.g. "~/ff/dreiding_den.par")

	suffix		A string which represents the job suffix. e.g. "tetramer"

	out_file	A string which represents the name of the output file.

	monomers	An integer which represents the number of PEIs that participates in a crosslinker.
			If # of bgf_file < monomers, then additional PEI structure files will be added by looking up files from 
			"/home/noische/research/dendrimer/structure/600/minEnergyOrder". You can change this directory by changing
			structure_dir variables in the script.

	Rg		A float which contains the distance between monomers.

	directory	A string which contains the pool of BGF files.

	silent		True/False flag for determining printing output messages.
	"""

	if monomers == 4:
		pass;
	elif monomers == 10:
		pass;
	elif monomers == 16:
		pass;
	elif monomers == 27:
		pass;
	elif monomers == 2:
		pass;
	else:
		nu.die("Sorry, %s crosslinking is not yet provided." % str(monomers))

	# variables
	chosen_pei = []
	extension = ".bgf"
	log_file = suffix + ".log"
	f_log_file = open(log_file, 'w')

	f_log_file.write("fourPEI.py: version " + str(version) + "\n")
	f_log_file.write("" + "\n")
	f_log_file.write("Job started at " + time.asctime(time.gmtime()) + " on " + os.environ["HOSTNAME"] + " by " + os.environ["USER"] + "\n")
	f_log_file.write("Command executed at " + os.getcwd() + "\n")
	f_log_file.write("Requested options: " + str(sys.argv) + "\n")
	f_log_file.write("" + "\n")

	# CONSTANTS
	#Rg = 6.0055	# average radius of gyration of PEI 600
	#Rg = 8.4088	# don't know?
	#Rg = 12.0100	# just feeling: actually Rg * 2
	f_log_file.write("Using the PEI intermolecular distance " + str(Rg) + "\n")

	# file preparation
	#structure_dir = "/home/noische/research/dendrimer/structure/600/minEnergyOrder"		# basic bgf repository
	structure_dir = directory
	curr_dir = os.path.abspath(".")
	temp_dir = os.path.join(curr_dir, suffix)
	pei_file = glob.glob(structure_dir + "/*.bgf")	# only extracts bgf filename.ext
	#pei_file = [os.path.basename(i) for i in glob.glob(structure_dir + "/*.bgf")]	# only extracts bgf filename.ext
	pei_name = [i.rstrip(".bgf\n")[0] for i in pei_file]				# only extracts bgf filename

	chosen_pei = bgf_file.split(" ")
	if chosen_pei == ['']: chosen_pei = []

	if len(chosen_pei) != monomers:
		for file in chosen_pei:
			if file in pei_file: pei_file.remove(file)	# prevents repeated selections of monomer
		for i in range(0, monomers - len(chosen_pei)):
			a = random.choice(pei_file)
			chosen_pei.append(a)
			pei_file.remove(a)	# prevents repeat. Now chosen_pei contains the files
	# REMARK: chosen_pei = ['a.bgf', 'b.bgf', 'c.bgf', '/home/noische/research/dendrimer/structure/600/minEnergyOrder/600_E025.bgf']

	# scratch preparation
	if not os.path.isdir(temp_dir):
		os.makedirs(temp_dir)
	else:
		#shutil.rmtree(temp_dir)
		os.makedirs(temp_dir)

	# copy files to the scratch
	for file in chosen_pei:
		shutil.copy(file, temp_dir)

	# move to the scratch directory
	os.chdir(temp_dir)

	chosen_pei = [os.path.basename(i) for i in chosen_pei]			# REMARK: ['a.bgf', 'b.bgf', 'c.bgf', '600_E018.bgf']
	chosen_pei_name = [i.partition(".bgf")[0] for i in chosen_pei]		# REMARK: ['a', 'b', 'c', '600_E018']

	f_log_file.write("Structure file directory: " + structure_dir + "\n")
	f_log_file.write("Current working directory: " + curr_dir + "\n")
	f_log_file.write("Temp file directory: " + temp_dir + "\n")
	f_log_file.write("Chosen PEI molecules: " + str(chosen_pei) + "\n")
	if not silent: print(chosen_pei)

	# open BGF structures as BgfFile
	myBGF = [ bgf.BgfFile(i) for i in chosen_pei ]

	# center BGF
	for structure in myBGF:
		structure = centerBGF.centerbgf(structure, 0, ff_file, "com_center", True)
	
	# move BGF for tetrahedron
	#moveBGF(myBGF[1], 0, Rg, Rg)
	#moveBGF(myBGF[2], Rg, 0, Rg)
	#moveBGF(myBGF[3], Rg, Rg, 0)

	# move BGF for planar
	if monomers == 4:
		# Tetrahedron
		bgf.moveBGF(myBGF[1], Rg, Rg, 0)
		bgf.moveBGF(myBGF[2], Rg, 0, Rg)
		bgf.moveBGF(myBGF[3], 0, Rg, Rg)
	elif monomers == 10:
		# Decahedron
		bgf.moveBGF(myBGF[1], Rg, 0, 0)
		bgf.moveBGF(myBGF[2], 2 * Rg, 0, 0)
		bgf.moveBGF(myBGF[3], 0.5 * Rg, 0.86603 * Rg, 0)
		bgf.moveBGF(myBGF[4], 1.5 * Rg, 0.86603 * Rg, 0)
		bgf.moveBGF(myBGF[5], Rg, 1.73205 * Rg, 0)
		bgf.moveBGF(myBGF[6], 0.5 * Rg, 0.28868 * Rg, 0.81650 * Rg)
		bgf.moveBGF(myBGF[7], 1.5 * Rg, 0.28868 * Rg, 0.81650 * Rg)
		bgf.moveBGF(myBGF[8], Rg, 1.15470 * Rg, 0.81650 * Rg)
		bgf.moveBGF(myBGF[9], Rg, 0.57735 * Rg, 1.63300 * Rg)
	elif monomers == 16:
		"""
		# < planar >
		bgf.moveBGF(myBGF[1], 0, Rg, 0)
		bgf.moveBGF(myBGF[2], 0, 2*Rg, 0)
		bgf.moveBGF(myBGF[3], 0, 3*Rg, 0)
		bgf.moveBGF(myBGF[4], Rg, 0, 0)
		bgf.moveBGF(myBGF[5], Rg, Rg, 0)
		bgf.moveBGF(myBGF[6], Rg, 2*Rg, 0)
		bgf.moveBGF(myBGF[7], Rg, 3*Rg, 0)
		bgf.moveBGF(myBGF[8], 2*Rg, 0, 0)
		bgf.moveBGF(myBGF[9], 2*Rg, Rg, 0)
		bgf.moveBGF(myBGF[10], 2*Rg, 2*Rg, 0)
		bgf.moveBGF(myBGF[11], 2*Rg, 3*Rg, 0)
		bgf.moveBGF(myBGF[12], 3*Rg, 0, 0)
		bgf.moveBGF(myBGF[13], 3*Rg, Rg, 0)
		bgf.moveBGF(myBGF[14], 3*Rg, 2*Rg, 0)
		bgf.moveBGF(myBGF[15], 3*Rg, 3*Rg, 0)
		"""
		# L shaped box
		bgf.moveBGF(myBGF[1], Rg, 2*Rg, 2*Rg)
		bgf.moveBGF(myBGF[2], 0, Rg, 0)
		bgf.moveBGF(myBGF[3], Rg, Rg, 0)
		bgf.moveBGF(myBGF[4], 0, 2*Rg, 0)
		bgf.moveBGF(myBGF[5], Rg, 2*Rg, 0)
		bgf.moveBGF(myBGF[6], Rg, 0, Rg)
		bgf.moveBGF(myBGF[7], Rg, Rg, Rg)
		bgf.moveBGF(myBGF[8], Rg, 2*Rg, Rg)
		bgf.moveBGF(myBGF[9], 0, 2*Rg, Rg)
		bgf.moveBGF(myBGF[10], 0, Rg, Rg)
		bgf.moveBGF(myBGF[11], 0, 0, Rg)
		bgf.moveBGF(myBGF[12], Rg, 0, 2*Rg)
		bgf.moveBGF(myBGF[13], Rg, Rg, 2*Rg)
		bgf.moveBGF(myBGF[14], 0, 0, 2*Rg)
		bgf.moveBGF(myBGF[15], 0, Rg, 2*Rg)
	elif monomers == 27:
		# cubic style crosslinking
		bgf.moveBGF(myBGF[1],   Rg,    0,    0)
		bgf.moveBGF(myBGF[2], 2*Rg,    0,    0)
		bgf.moveBGF(myBGF[3],    0,   Rg,    0)
		bgf.moveBGF(myBGF[4],   Rg,   Rg,    0)
		bgf.moveBGF(myBGF[5], 2*Rg,   Rg,    0)
		bgf.moveBGF(myBGF[6],    0, 2*Rg,    0)
		bgf.moveBGF(myBGF[7],   Rg, 2*Rg,    0)
		bgf.moveBGF(myBGF[8], 2*Rg, 2*Rg,    0)
		bgf.moveBGF(myBGF[9],    0,    0,   Rg)
		bgf.moveBGF(myBGF[10],   Rg,    0,   Rg)
		bgf.moveBGF(myBGF[11], 2*Rg,    0,   Rg)
		bgf.moveBGF(myBGF[12],    0,   Rg,   Rg)
		bgf.moveBGF(myBGF[13],   Rg,   Rg,   Rg)
		bgf.moveBGF(myBGF[14], 2*Rg,   Rg,   Rg)
		bgf.moveBGF(myBGF[15],    0, 2*Rg,   Rg)
		bgf.moveBGF(myBGF[16],   Rg, 2*Rg,   Rg)
		bgf.moveBGF(myBGF[17], 2*Rg, 2*Rg,   Rg)
		bgf.moveBGF(myBGF[18],    0,    0, 2*Rg)
		bgf.moveBGF(myBGF[19],   Rg,    0, 2*Rg)
		bgf.moveBGF(myBGF[20], 2*Rg,    0, 2*Rg)
		bgf.moveBGF(myBGF[21],    0,   Rg, 2*Rg)
		bgf.moveBGF(myBGF[22],   Rg,   Rg, 2*Rg)
		bgf.moveBGF(myBGF[23], 2*Rg,   Rg, 2*Rg)
		bgf.moveBGF(myBGF[24],    0, 2*Rg, 2*Rg)
		bgf.moveBGF(myBGF[25],   Rg, 2*Rg, 2*Rg)
		bgf.moveBGF(myBGF[26], 2*Rg, 2*Rg, 2*Rg)
	elif monomers == 2:
		bgf.moveBGF(myBGF[1], Rg, 0, 0)

	# change resnumber of monomers
	for index, structure in enumerate(myBGF):
		for atom in structure.a:
			#atom.rNo = index + 1
			#atom.rNo += (index + 1) * 10
			atom.rNo += (index + 1) * 100

	# merge BGF for any number of monomers
	"""
	myBGF[0] = myBGF[0].merge(myBGF[1], True)
	myBGF[0] = myBGF[0].merge(myBGF[2], True)
	myBGF[0] = myBGF[0].merge(myBGF[3], True)
	"""
	for i in range(1, monomers):
		myBGF[0] = myBGF[0].merge(myBGF[i], True)

	# record BGF properties
	myBGF[0].REMARK.insert(0, "Crosslinked structures: " + str(chosen_pei))
	myBGF[0].REMARK.insert(0, "Crosslinked by " + os.path.basename(sys.argv[0]) + " by " + os.environ["USER"] + " on " + time.asctime(time.gmtime()))

	# save BGF
	myBGF[0].saveBGF(suffix + ".bgf")

	# crosslink
	#if yes_crosslink:
	#	f_log_file.write("Entering crosslinking.\n")
	#	crosslink.crosslink(suffix + ".bgf", ff_file, suffix, nodes, criteria, t, out_file)
	#else:
	#	f_log_file.write("-M option activated. Ignore crosslinking process.\n")

	# file copy
	shutil.copy(os.path.join(temp_dir, suffix + ".bgf"), curr_dir)
	#os.system("cp -rp " + temp_dir + " " + curr_dir)

	f_log_file.close()

	### end of fourPEI.py

if __name__ == '__main__':

	option = ""; args = ""; bgf_file = ""; ff_file = ""; out_file = "";
	suffix = ""; monomers = ""; Rg = 0.0; directory = "";
	usage = """
Usage: fourPEI.py -h -b "bgfFiles" -f forcefield -s suffix -n nodes -c xlnk_dist_criteria -t n_cycles -m monomers -M yes/no

Options are:
	-h	Print help.
	-b	A series of BGF monomers in quotes("").
	-m	Number of monomers that will be crosslinked. 
		If the number of designated BGF file with -b option is less than the number,
		BGF files in the pool (/home/noische/research/dendrimer/structure/600/minEnergyOrder) will be automatically chosen.
		to meet the number in this option.
	-s	Suffix of the job script (common in this script and crosslink.py script)
	-o	Output file. (DEFAULT suffix_output.log)
	-m	Number of monomers: 2, 4, 10, 16, 27 are available.
	-d	Distance between PEI molecules. (DEFAULT 6.0055)
	-D	Directory which contains monomers in BGF.
	"""

	if len(sys.argv) < 2:
		print(usage);
		sys.exit(0)

	options, args = getopt.getopt(sys.argv[1:], 'hb:f:s:o:m:d:D:', ['help','bgf=','forcefield=','suffix=','output=','monomers=','distance=','directory='])
	for option, value in options:
		if option in ('-h', '--help'):
			print usage; sys.exit(0)
		elif option in ('-b', '--bgf'):
			bgf_file = value
		elif option in ('-f', '--forcefield'):
			ff_file = value
		elif option in ('-s', '--suffix'):
			suffix = value
		elif option in ('-o', '--output'):
			out_file = value
		elif option in ('-m', '--monomers'):
			monomers = int(value)
		elif option in ('-d', '--distance'):
			Rg = float(value);
		elif option in ('-D', '--directory'):
			directory = value
		elif option in (''):
			print(usage); sys.exit(0)

	# default settings
    if not suffix: suffix = "fourpei"
	if not out_file: out_file = suffix + "_output.log"
	if not ff_file: ff_file = "/home/noische/ff/DREIDING2.21.ff"
	if not Rg: Rg = 6.0055		# Average Radius of Gyration of PEI 600 monomers
	if not directory: directory = "/home/noische/research/dendrimer/structure/600/minEnergyOrder"		# basic bgf repository
	fourPEI(bgf_file, ff_file, suffix, out_file, monomers, Rg, directory, silent=False)
