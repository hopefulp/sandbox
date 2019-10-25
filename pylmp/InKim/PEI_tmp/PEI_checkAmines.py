#!/home/noische/python
"""
test_isomorphism.py
Original: Dec 28 2011 In Kim
"""

# Python Modules
import sys
import os
import string
import random
import time
import getopt
import glob

# Custom Modules
sys.path.append("/qcfs/noische/scripts")
import bgf
import bgftools
import nutils as nu

# Globals
# 120102
# 140627 added time estimation
#	 modified as a module to use in generatePolymers.py

version = '140627'
usage = 'PEI_checkAmines.py bgf_directory'

def check_amines(directory):
	"""
	check isomorphism of 'bgf_file' in all files in 'directory'
	"""

	#initialize
	structure_dir = os.path.abspath(directory)
	pei_file = glob.glob(structure_dir + "/*.bgf")
	pei_file.sort()
	n_pei_file = len(pei_file)
	print("The script will check " + str(n_pei_file) + " files in the directory " + structure_dir)

	joblist = [];	# contains file pairlists
	n_joblist = 0;	# number of total jobs
	count = 0;	# job counter
	t1 = time.time(); t2 = 0;	# time progess


	for job in pei_file:
		# initialize
		t2 = time.time()
		elapsed = t2 - t1
		estimated = elapsed * (n_joblist - count)
	
		# open BGF
		myBGF = bgf.BgfFile(job)

		print(bgftools.getAmineGroupInfo(myBGF))
	
		# display the process
		#sys.stdout.write("\rProgress: " + "{0:>8d}".format(count) + " / " + str(n_joblist) + " (" + str(estimated) + " sec left)"); sys.stdout.flush()

	# completing
	return True;

	### end of check_amines()

if len(sys.argv) < 2:
	print(usage)
	sys.exit(0)

check_amines(sys.argv[1])
