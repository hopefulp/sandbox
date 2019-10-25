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
version = '150921 kdft'

def do_work(bgf_file, out_file, simple, silent=True):

	#initialize
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

    # get Wiener Index
    index = 0;
    for key in d_dist.iterkeys():
        for key2 in d_dist[key].iterkeys():
            index += d_dist[key][key2];
    index = index / 2.0

    return index

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
