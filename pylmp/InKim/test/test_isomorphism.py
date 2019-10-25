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
f_out_file = open("isomorphism_check.dat", "w")

def test_isomorphism(bgf_file1, bgf_file2, simple, silent=True):
	"""
def test_isomorphism():

Function Parameters:
	bgf_file	A string of filename or BgfFile class.

	"""
	# initialize
	pei1 = 0; pei2 = 0; myBGF1 = 0; myBGF2 = 0;

	# open BGF
	if isinstance(bgf_file1, bgf.BgfFile):
		myBGF1 = bgf_file1
	else:
		if not silent: print("opening bgf file 1.. " + str(bgf_file1))
		myBGF1 = bgf.BgfFile(bgf_file1)

	if isinstance(bgf_file2, bgf.BgfFile):
		myBGF2 = bgf_file2
	else:
		if not silent: print("opening bgf file 2.. " + str(bgf_file2))
		myBGF2 = bgf.BgfFile(bgf_file2)

	if simple:
		# remove hydrogens only in carbon atoms
		remove_aNo1 = []; remove_aNo2 = [];
		for atom in myBGF1.a:
			if "C_" in atom.ffType:
				n_H = [];
				for aNo in atom.CONECT:
					if myBGF1.getAtom(aNo).is_hydrogen(): n_H.append(aNo)
				print("aNo: %s, n_H: %s in %s (1)" % (atom.aNo, n_H, bgf_file1))
			if len(n_H) != 3:
				myBGF1.delAtoms(n_H)
				myBGF1.renumber()
		#print("saving mybgf1"); myBGF1.saveBGF(bgf_file1 + "mod")
		pei1 = myBGF1

		for atom in myBGF2.a:
			if "C_" in atom.ffType:
				n_H = [];
				for aNo in atom.CONECT:
					if myBGF2.getAtom(aNo).is_hydrogen(): n_H.append(aNo)
				print("aNo: %s, n_H: %s in %s (2)" % (atom.aNo, n_H, bgf_file2))
			if len(n_H) != 3:
				myBGF2.delAtoms(n_H)
				myBGF2.renumber()
		#print("saving mybgf2"); myBGF2.saveBGF(bgf_file2 + "mod")
		pei2 = myBGF2
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


def do_work(qin, qout, simple):
	while True:
		x = qin.get()
		if x == None:
			break;
		elif qin.qsize() < 2:
			break;
		else:
			# initialize
			bgf_file1 = x[0]; bgf_file2 = x[1];
			pei1 = 0; pei2 = 0; myBGF1 = 0; myBGF2 = 0; n_H = [];
			pid = os.getpid()	# process id
		
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
		
			# to check reduced structures
			#pei1.saveBGF(bgf_file1 + "mod")
			#pei2.saveBGF(bgf_file2 + "mod")
	
			# convert connection into dictionary
			d_graph1 = bgftools.getConnectionDict(pei1)
			d_graph2 = bgftools.getConnectionDict(pei2)
		
			# convert dictionary into Graph
			G1 = nx.Graph(d_graph1)
			G2 = nx.Graph(d_graph2)
		
			# check isomorphism
			GM = isomorphism.GraphMatcher(G1, G2)
			result = GM.is_isomorphic()
	
			# collect outputs
			qout.put((pid, os.path.basename(bgf_file1), os.path.basename(bgf_file2), result))
	
			# display the process
			sys.stdout.write("\rNumber of leftover: " + "{0:>8d}".format(qin.qsize()) + ", Number in qout: " + "{0:>8d}".format(qout.qsize())); sys.stdout.flush()

	return 1;


def main(directory, nprocs, out_file, simple, silent):

	#initialize
	structure_dir = os.path.abspath(directory)
	curr_dir = os.path.abspath(".")
	pei_file = glob.glob(structure_dir + "/*.bgf")
	pei_file.sort()
	n_pei_file = len(pei_file)

	qin = multiprocessing.Queue()	# Queue for two bgf files
	qout = multiprocessing.Queue()	# Queue for output

	for i in range(0, n_pei_file):
		for j in range(i, n_pei_file):
			if not i == j:
				qin.put([pei_file[i], pei_file[j]]);

	if not silent: print("Files on Queue. Processing.. Total number of calculations on Queue: " + "{0:>8d}".format(qin.qsize()) )

	workers = [];
	for i in range(nprocs):
		qin.put(None)	# sentinel

	for i in range(nprocs):
		tmp = multiprocessing.Process(target=do_work, args=(qin, qout, simple, ))
		tmp.start()
		workers.append(tmp)

	try:
		for w in workers:
			w.join();

	except KeyboardInterrupt:
		print("Keyboard Break.. Quitting Job.")
		for w in workers:
			w.terminate()
		qout.put(None)	# sentinel
		writeQueue(qout, out_file)
		return False;

	except:
		print("Unknown Error Occurred.. Quitting Job.")
		for w in workers:
			w.terminate()
		qout.put(None)	# sentinel
		writeQueue(qout, out_file)
		return False;
		

	# record the result
	qout.put(None)	# sentinel
	writeQueue(qout, out_file)
	return True;

	### end of main()


def writeQueue(queue, file):
	f = open(file, 'w')
	while True:
		try:
			x = queue.get(block = False)

			# output: (pid, os.path.basename(bgf_file1), os.path.basename(bgf_file2), result)
			output = "(pid %s)\t%s\t%s\t%s\n" % x
			f.write(output)
		except:
			break;
	f.close()

	### end of writeQueue


if __name__ == '__main__':
	multiprocessing.freeze_support()

	option = ""; args = ""; bgf_file1 = ""; bgf_file2 = ""; directory = ""; out_file = ""; nprocs = 0; silent = False; simple = False;
	usage = """
Usage: test_isomorphism.py -b "bgfFiles" -s size -n monomers -o output

Options are:
	-d	Input BGF directory
	-r	Activate reduced mode
	-n	Number of processes calculation
	-s	Silent mode (no printout)
	"""

	if len(sys.argv) < 2:
		print(usage); sys.exit(0)

	options, args = getopt.getopt(sys.argv[1:], 'hb:c:d:rn:o:s', ['help','1=','2=','directory=','simple=','nprocs=','output=','silent='])
	for option, value in options:
		if option in ('-h', '--help'):
			print usage; sys.exit(0)
		elif option in ('-b', '--1'):
			bgf_file1 = value
		elif option in ('-c', '--2'):
			bgf_file2 = value
		elif option in ('-d', '--directory'):
			directory = value
		elif option in ('-r', '--reduce'):
			simple = True
		elif option in ('-n', '--nprocs'):
			nprocs = int(value)
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
	a = main(directory, nprocs, out_file, simple, silent)

	if a:
		print("\nJob Completed")
	else:
		print("\nError")
