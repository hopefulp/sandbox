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
sys.path.append("/home/noische/scripts")
sys.path.append("/home/noische/script")
import bgf
import bgftools
import nutils as nu
import networkx as nx
from networkx.algorithms import isomorphism

# Globals
# 120102
# 140627 added time estimation
#	 modified as a module to use in generatePolymers.py

version = '140627'
usage = """
Check isomorphism of 'bgf_file' in all files in 'directory' in simple mode.
usage: check_isomorphism.py directory bgf_file
"""

def check_isomorphism(directory, bgf_file, simple):
	"""
	check isomorphism of 'bgf_file' in all files in 'directory'
	"""

	#initialize
	structure_dir = os.path.abspath(directory)
	curr_dir = os.path.abspath(".")
	pei_file = glob.glob(structure_dir + "/*.bgf")
	for i in pei_file:
		if bgf_file in i:
			pei_file.remove(i)	# remove self bgf_file from pei_file

	pei_file.sort()
	n_pei_file = len(pei_file)
	print("The script will compare " + str(n_pei_file) + " files in the directory " + curr_dir)

	joblist = [];	# contains file pairlists
	n_joblist = 0;	# number of total jobs
	count = 0;	# job counter
	t1 = time.time(); t2 = 0;	# time progess

	print("Queueing Jobs..")
	print(n_pei_file)
	if n_pei_file == 0:
		return True;

	for i in range(0, n_pei_file):
		joblist.append([bgf_file, pei_file[i]]);
	n_joblist = len(joblist)

	myBGF1 = bgf.BgfFile(bgf_file)
	if simple:
		# remove hydrogens only in carbon atoms
		pei1 = bgftools.getBackbone(bgf_file1)
	else:
		pei1 = bgf_file1
	

	for job in joblist:
		# initialize
		bgf_file2 = job[1];		# filenames
		pei1 = 0; pei2 = 0; myBGF2 = 0; n_H = [];
		t2 = time.time()
		elapsed = t2 - t1
		estimated = elapsed * (n_joblist - count)
	
		# open BGF
		myBGF2 = bgf.BgfFile(bgf_file2)
	
		if simple:
			# remove hydrogens only in carbon atoms
			pei2 = bgftools.getBackbone(bgf_file2)
		else:
			# toss whole structures
			pei2 = bgf_file2
	
		# convert connection into dictionary
		d_graph1 = bgftools.getConnectionDict(pei1)
		d_graph2 = bgftools.getConnectionDict(pei2)
	
		# convert dictionary into Graph
		G1 = nx.Graph(d_graph1)
		G2 = nx.Graph(d_graph2)
	
		# check isomorphism
		GM = isomorphism.GraphMatcher(G1, G2)
		result = GM.is_isomorphic()

		if not result:
			return bgf_file2;
			#print bgf_file2

		# count
		count += 1;
	
		# display the process
		#sys.stdout.write("\rProgress: " + "{0:>8d}".format(count) + " / " + str(n_joblist) + " (" + str(estimated) + " sec left)"); sys.stdout.flush()

	# completing
	return True;

	### end of check_isomorphism()


if __name__ == '__main__':
	if len(sys.argv) < 2:
		print(usage)
		sys.exit(0)

	check_isomorphism(sys.argv[1], sys.argv[2], True)
