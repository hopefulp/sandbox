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

version = '140627'

def do_work(directory, out_file, simple, silent=True):

	#initialize
	structure_dir = os.path.abspath(directory)
	curr_dir = os.path.abspath(".")
	pei_file = glob.glob(structure_dir + "/*.bgf")
	pei_file.sort()
	n_pei_file = len(pei_file)
	f_out_file = open(out_file, 'w')	# file for result 
	f_out_file.write(str(sys.argv))
	joblist = [];	# contains file pairlists
	n_joblist = 0;	# number of total jobs
	count = 0;	# job counter
	t1 = time.time(); t2 = 0;	# time progess

	if not silent: print("Queueing Jobs..")
	for i in range(0, n_pei_file):
		for j in range(i, n_pei_file):
			if not i == j:
				joblist.append([pei_file[i], pei_file[j]]);
	n_joblist = len(joblist)

	for job in joblist:
		# initialize
		bgf_file1 = job[0]; bgf_file2 = job[1];		# filenames
		pei1 = 0; pei2 = 0; myBGF1 = 0; myBGF2 = 0; n_H = [];
		t2 = time.time()
		elapsed = t2 - t1
		estimated = elapsed * (n_joblist - count)
	
		# open BGF
		myBGF1 = bgf.BgfFile(bgf_file1)
		myBGF2 = bgf.BgfFile(bgf_file2)
	
		if simple:
			# remove hydrogens only in carbon atoms
			remove_aNo1 = []; remove_aNo2 = [];
			for atom in myBGF1.a:
				if "C_" in atom.ffType:
					n_H = [];
					for aNo in atom.CONECT:
						if myBGF1.getAtom(aNo).is_hydrogen(): n_H.append(myBGF1.a2i[aNo])
					if len(n_H) != 3:
						n_H.sort(); n_H.reverse(); myBGF1.delAtoms(n_H); myBGF1.renumber()
	
			for atom in myBGF2.a:
				if "C_" in atom.ffType:
					n_H = [];
					for aNo in atom.CONECT:
						if myBGF2.getAtom(aNo).is_hydrogen(): n_H.append(myBGF2.a2i[aNo])
					if len(n_H) != 3:
						n_H.sort(); n_H.reverse(); myBGF2.delAtoms(n_H); myBGF2.renumber()
	
			# toss reduced BGFs to networkx
			pei1 = myBGF1; pei2 = myBGF2
	
		else:
			# toss whole structures
			pei1 = bgf_file1
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

		# count
		count += 1;
	
		# write outputs
		output = os.path.basename(bgf_file1) + "\t" + os.path.basename(bgf_file2) + "\t" + str(result) + "\n"
		f_out_file.write(output)
	
		# display the process
		if not silent: sys.stdout.write("\rProgress: " + "{0:>8d}".format(count) + " / " + str(n_joblist) + " (" + str(estimated) + " sec left)"); sys.stdout.flush()

	# completing
	f_out_file.close()

	return 1;


if __name__ == '__main__':

	option = ""; args = ""; bgf_file1 = ""; bgf_file2 = ""; directory = ""; out_file = ""; nprocs = 0; silent = False; simple = False;
	usage = """
Usage: test_isomorphism.py -b "bgfFiles" -s size -n monomers -o output

Options are:
	-d	Input BGF directory
	-r	Activate reduced mode
	-s	Silent mode (no printout)
	"""

	if len(sys.argv) < 2:
		print(usage); sys.exit(0)

	options, args = getopt.getopt(sys.argv[1:], 'hd:ro:s', ['help','directory=','simple=','output=','silent='])
	for option, value in options:
		if option in ('-h', '--help'):
			print usage; sys.exit(0)
		elif option in ('-d', '--directory'):
			directory = value
		elif option in ('-r', '--reduce'):
			simple = True
		elif option in ('-o', '--output'):
			out_file = value
		elif option in ('-s', '--silent'):
			silent = True
		elif option in (''):
			print(usage); sys.exit(0)

	# default settings
	if out_file == "":
		out_file = "isomorphism.%s.dat" % time.strftime("%Y%m%d-%I%M", time.localtime())

	#test_isomorphism(bgf_file1, bgf_file2, silen=False)
	a = do_work(directory, out_file, simple, silent)

	if a:
		print("\nJob Completed")
	else:
		print("\nError")
	
