#!/home/noische/python
#!/opt/applic/epd/bin/python
"""
PEI_calcWienerIndex.py
Original: Jul 11 2012 In Kim
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
import bgf
import bgftools
import nutils as nu
import networkx as nx

# Globals
version = '120102 kdft'

def do_work(directory, out_file, simple, silent=True):

	#initialize
	structure_dir = os.path.abspath(directory)
	curr_dir = os.path.abspath(".")
	pei_file = glob.glob(structure_dir + "/*.bgf")
	pei_file.sort()
	n_pei_file = len(pei_file)
	f_out_file = open(out_file, 'w')	# file for result 
	joblist = [];	# contains filelists
	n_joblist = 0;	# number of total jobs
	count = 0;	# job counter
	l_index = [];	# for index

	f_out_file.write(str(sys.argv[0]) + " version " + str(version) + "\n")
	f_out_file.write("" + "\n")
	f_out_file.write("Job started at " + time.asctime(time.gmtime()) + " on " + os.environ["HOSTNAME"] + " by " + os.environ["USER"] + "\n")
	f_out_file.write("Command executed at " + os.getcwd() + "\n")
	f_out_file.write("Requested options: " + str(sys.argv) + "\n")
	f_out_file.write("" + "\n")

	if not silent: print("Queueing Jobs..")
	for i in range(0, n_pei_file):
		joblist.append(pei_file[i]);
	n_joblist = len(joblist)

	for job in joblist:
		# initialize
		bgf_file = job;
		myPEI = bgf.BgfFile(); n_H = [];
	
		# open BGF
		myBGF = bgf.BgfFile(bgf_file)
	
		# remove hydrogens 
		pei = bgftools.getBackbone(myBGF, 0)

		# convert connection into dictionary
		d_graph = bgftools.getConnectionDict(pei)
	
		# convert dictionary into Graph
		G = nx.Graph(d_graph)
	
		# calculate the all shortest length path
		#d_dist = nx.floyd_warshall(G)	# seems a problem
		d_dist = nx.all_pairs_shortest_path_length(G)

		# count
		count += 1;

		# get Wiener Index
		index = 0;
		for key in d_dist.iterkeys():
			for key2 in d_dist[key].iterkeys():
				index += d_dist[key][key2];
		index = index / 2.0

		# find duplicates
		if index in l_index:
			output = os.path.basename(bgf_file) + "\t" + str(index) + "\t" + "duplicate" + "\n"
		else:
			output = os.path.basename(bgf_file) + "\t" + str(index) + "\n"
	
		# for index checking
		l_index.append(index)

		# write outputs
		f_out_file.write(output)
	
		# display the process
		if not silent: sys.stdout.write("\rProgress: " + "{0:>8d}".format(count) + " / " + str(n_joblist)); sys.stdout.flush()

	# check duplicates
	tmp=[]; dup=[];
	for i in l_index:
		if i in tmp:
			dup.append(i)
		else:
			tmp.append(i)

	f_out_file.write("\nDuplicated Wiener Indices:\n")
	f_out_file.write(str(dup) + "\n")
	f_out_file.close()

	return 1;

	### end of function


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

	options, args = getopt.getopt(sys.argv[1:], 'hd:o:s', ['help','directory=','output=','silent'])
	for option, value in options:
		if option in ('-h', '--help'):
			print usage; sys.exit(0)
		elif option in ('-d', '--directory'):
			directory = value
		elif option in ('-o', '--output'):
			out_file = value
		elif option in ('-s', '--silent'):
			silent = True
		elif option in (''):
			print(usage); sys.exit(0)

	# default settings
	if out_file == "":
		out_file = "WienerIndex.%s.dat" % time.strftime("%Y%m%d-%I%M", time.localtime())

	#test_isomorphism(bgf_file1, bgf_file2, silen=False)
	a = do_work(directory, out_file, silent)

	if a:
		print("\nJob Completed")
	else:
		print("\nError")
