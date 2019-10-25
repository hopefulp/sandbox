#!/home/noische/program/python27/bin/python
"""
poresize.py
Original: Apr 26 2011 In Kim
"""

# Python Modules
import sys
import os
import string
import random
import time
import getopt
import numpy as np
import scipy.io

# Custom Modules
import bgf
import nutils as nu

# Globals
version = '110501'

def poresize(bgf_file, out_file, grid = 0.25, thresh = 0.25, silent=False):
	"""
def poresize():
	Calculates pore sizes in the BGF molecule

Function Parameters:
	bgf_file	A string which contains a PEI information in a BGF format. 
			(e.g. "file1.bgf file2.bgf file3.bgf ...")
	"""

	# initialize
	xyz = [];
	rvdw = dict();
	ubound = 10; lbound = 0;

	# open bgf
	myBGF = bgf.BgfFile(bgf_file)

	# read van der waals radii
	rvdw['H'] = 1.2
	rvdw['C'] = 1.7
	rvdw['N'] = 1.55
	rvdw['O'] = 1.52
	
	# copy positions: [ x y z vdw ]
	for atom in myBGF.a:
		atomtype = atom.ffType[:-1].strip("_")
		xyz.append([ float(atom.x), float(atom.y), float(atom.z), rvdw[atomtype] ])

	xyz = np.array(xyz)

	# get bgf size
	size = np.array(bgf.getBGFSize(myBGF))
	np.around(size, decimals=1)

	# make grid
	grid_x = np.arange(size[0], size[1], grid)
	grid_y = np.arange(size[2], size[3], grid)
	grid_z = np.arange(size[4], size[5], grid)
	len_grid_x = len(grid_x)
	len_grid_y = len(grid_y)
	len_grid_z = len(grid_z)
	poresize = np.zeros((len_grid_x,len_grid_y,len_grid_z))

	ix = 0; iy = 0; iz = 0; t1 = 0; t2 = 0;

	# for each grid
	for ix, x in enumerate(grid_x):
		t1 = time.time();

		for iy, y in enumerate(grid_y):

			for iz, z in enumerate(grid_z):

				ubound = 10; lbound = 0; minubound = 10;

				while ubound - lbound > thresh:
				
					test_r = random.uniform(lbound, ubound);

					flag = True;
					for atom in xyz:
						value = (x - atom[0])**2 + (y - atom[1])**2 + (z - atom[2])**2 - atom[3]**2 - test_r**2
						if value < 0: 
							flag = False
							if minubound > test_r: minubound = test_r

					if flag:
						lbound = test_r
					else:
						ubound = minubound

					poresize[ix, iy, iz] = test_r

		t2 = time.time();
		print(str(ix) + " out of " + str(len_grid_x) + "( elapsed time for the step: " + str(t2 - t1) + " sec)")

	# store them as a matrix
	scipy.io.savemat(out_file, {'m':poresize})
	
	# end

if __name__ == '__main__':

	option = ""; args = ""; bgf_file = ""; grid = ""; out_file = "";
	suffix = ""; nodes = ""; criteria = 0.0; t = 0; monomers = "";
	usage = """
Usage: poresize.py -b "bgfFiles" -f forcefield -s suffix -n nodes -c xlnk_dist_criteria -t n_cycles -m monomers

Options are:
	-b	A series of BGF monomers in quotes("").
	-o	Output file. (DEFAULT suffix_output.log)
	"""

	options, args = getopt.getopt(sys.argv[1:], 'hb:f:s:n:c:t:o:m:', ['help','bgf=','forcefield=','suffix=','nodes=','criteria=','time=','output=','monomers='])
	for option, value in options:
		if option in ('-h', '--help'):
			print usage; sys.exit(0)
		elif option in ('-b', '--bgf'):
			bgf_file = value
		elif option in ('-o', '--output'):
			out_file = value
		elif option in (''):
			print usage; sys.exit(0)

	# default settings
	if not out_file: out_file = bgf_file[:-4] + "_poresize.mat"

	poresize(bgf_file, out_file)

