#!/opt/applic/epd/bin/python
#/home/noische/program/python27/bin/python
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
#from multiprocessing import Process, Queue
import multiprocessing
import Queue
import errno

# Custom Modules
sys.path.append("/home/noische/scripts")
sys.path.append("/home/noische/script")
import bgf
import bgftools
import nutils as nu
import networkx as nx
from networkx.algorithms import isomorphism

# Globals
version = '120102'
simple = False;
f_out_file = open("isomorphism_check.dat", "w")
output = [];

def test_isomorphism(bgf_file1, bgf_file2, simple, silent=True):
	"""
def test_isomorphism():

Function Parameters:
	bgf_file	A string of filename or BgfFile class.

	"""
	# initialize
	pei1 = 0; pei2 = 0;

	# open BGF
	if isinstance(bgf_file1, bgf.BgfFile):
		myBGF = bgf_file1
	else:
		if not silent: print("opening bgf file 1.. " + str(bgf_file1))
		myBGF = bgf.BgfFile(bgf_file1)

	if isinstance(bgf_file2, bgf.BgfFile):
		myBGF = bgf_file2
	else:
		if not silent: print("opening bgf file 2.. " + str(bgf_file2))
		myBGF = bgf.BgfFile(bgf_file2)

	if simple:
		# reduce hydrogens
		pei1 = bgftools.getBackbone(bgf_file1)
		pei2 = bgftools.getBackbone(bgf_file2)
	else:
		pei1 = bgf_file1
		pei2 = bgf_file2

	# convert connection into dictionary
	d_graph1 = bgftools.getConnectionDict(pei1)
	d_graph2 = bgftools.getConnectionDict(pei2)
	if not silent: print("structures are converted to graph.")

	G1 = nx.Graph(d_graph1)
	G2 = nx.Graph(d_graph2)
	if not silent: print("graphs are loaded.")

	# check isomorphism
	if not silent: print("now comparing two graphs..")
	GM = isomorphism.GraphMatcher(G1, G2)
	result = GM.is_isomorphic()

	if not silent:
		print(result)
	else:
		return result;

	### end of test_isomorphism


def do_work(l, q):
	while True:
		x = q.get()
		if x == "" or q.qsize == 0:
			break;
		else:
			#result = test_isomorphism(x[0], x[1], simple)
			#l.acquire()
			#f_out_file.write(str(os.path.basename(x[0])) + "\t" + str(os.path.basename(x[1])) + "\t" + str(result) + "\n")
			#print(str(os.path.basename(x[0])) + "\t" + str(os.path.basename(x[1])) + "\t" + str(result))
			#l.release()
			sys.stdout.write("\rNumber of leftover: " + "{0:8d}".format(q.qsize())); sys.stdout.flush()

def main(directory, simple):

	#initialize

	structure_dir = os.path.abspath(directory)
	curr_dir = os.path.abspath(".")
	pei_file = glob.glob(structure_dir + "/*.bgf")
	n_pei_file = len(pei_file)

	q = multiprocessing.Queue()

	for i in range(0, n_pei_file):
		for j in range(i, n_pei_file):
			if not i == j:
				q.put([pei_file[i], pei_file[j]]);
	q.put("")	# sentinel

	#f_out_file.write(str(structure_dir) + "\n")
	print("Files on Queue. Processing..")

	l = multiprocessing.Lock()
	processes = [multiprocessing.Process(target=do_work, args=(l,q,)) for i in range(4)]

	for p in processes:
		p.start()

	notintr = False;
	for p in processes:
		"""
		# from here
		while not notintr:
			try:
				p.join()
				notintr = True;
			except OSError, ose:
				if ose.errno != errno.EINTR:
					raise ose
		# to here
		"""
		p.join()
	"""
	workers = [];
	for i in range(4):
		tmp = multiprocessing.Process(target=do_work, args=(q,))
		tmp.start()
		workers.append(tmp)

	for w in workers:
		w.join()
	"""

	for i in output:
		f_out_file.write(i)

	sys.exit(0)

if __name__ == '__main__':
	multiprocessing.freeze_support()

	option = ""; args = ""; bgf_file1 = ""; bgf_file2 = ""; directory = ""; 
	usage = """
Usage: test_isomorphism.py -b "bgfFiles" -s size -n monomers -o output

Options are:
	-d	Input BGF directory.
	-s	Activate simple mode.
	"""

	if len(sys.argv) < 2:
		print(usage); sys.exit(0)

	options, args = getopt.getopt(sys.argv[1:], 'hb:c:d:s', ['help','1=','2=','directory=','simple='])
	for option, value in options:
		if option in ('-h', '--help'):
			print usage; sys.exit(0)
		elif option in ('-b', '--1'):
			bgf_file1 = value
		elif option in ('-c', '--2'):
			bgf_file2 = value
		elif option in ('-d', '--directory'):
			directory = value
		elif option in ('-s', '--simple'):
			simple = True
		elif option in (''):
			print(usage); sys.exit(0)

	# default settings

	#test_isomorphism(bgf_file1, bgf_file2, silen=False)
	main(directory, simple)

