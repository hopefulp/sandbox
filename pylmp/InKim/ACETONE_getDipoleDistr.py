#!/home/noische/program/python27/bin/python

import sys, re, string, getopt, optparse, math, time, pprint 
from os import popen
import bgf
import bgftools
import numpy as np
import operator
import copy
import nutils as nu
import datatools
import dreiding
import cPickle as pickle
import pprint

import dump

option = ""; args = ""; bgf_file = ""; trj_file = ""; out_file = ""; 
version = "130221"

#-----------------
# calculate dipole moments from the LAMMPS trajectory
# http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_moment.py
# last update: 2013/02/06
#_________________
def get_line(trj_file):
	"""
	Open a lazy file stream. This is considered to increase a file reading speed.
	"""
	with open(trj_file) as file:
		for i in file:
			yield i


def densityProfile(bgf_file, trj_file, ff_file, pickle_file, out_file, avg_timestep, silent=False):

	nu.warn("LAMMPS trajectory with NPT simulations will give you the wrong result.")

	### init
	timestep = 0; l_timestep = []; line = []; n_header = 0;
	t1 = 0; t2 = 0; # clock
	atominfo = dict();	# atom data extracted from ff_file
	result = dict();
	axis = 2;		# 1: x-axis, 2: y-axis, 3: z-axis


	### open files
	myBGF = bgf.BgfFile(bgf_file)
	myTRJ = open(trj_file)
	myTRJ.seek(0)
	myPickle = open(pickle_file)
	f_debug = open("debug.dat", 'w')
	if not silent: print("Pickling..")
	Dipole = pickle.load(myPickle)	# dipole data


	### read residues from bgf_file
	residue = set();	# kind of residues in BGF file
	dict_residue = dict();	# stores residue numbers per each residue. RES1: [1, 2, ..], RES2: [6, 7, ..]
	for i in myBGF.a:
		rname = string.strip(i.rName)
		residue.add(rname)
	output = "";
	for i in residue:
		output += i + " ";
		dict_residue[i] = [];

	if not silent: print("Found " + str(len(residue)) + " residues in BGF file: " + str(output))
	residue = list(residue)
	#residue.append('TOTAL')

	n_residue = len(residue)	# number of residues (including total)


	### bookkeep residue numbers for each residue
	molecules = bgftools.getMoleculeList(myBGF)
	dict_rNo2rName = dict();

	for molecule in molecules:
		# get an atom
		atom = myBGF.getAtom(molecule[0])
		atomresname = string.strip(atom.rName)
		atomresno = atom.rNo

		# check if all molecule has same residue name and numbers
		for ano in molecule:
			atom2 = myBGF.getAtom(ano)
			temp_rno = atom2.rNo
			temp_rname = string.strip(atom2.rName)
			if temp_rno != atomresno or temp_rname != atomresname:
				nu.die("Different residue name or residue number in a same molecule: " + str(atom2.aNo))

		# record resid for resnames
		dict_residue[atomresname].append(atom.rNo)
		dict_rNo2rName[atom.rNo] = atomresname

	#print(dict_residue)
	#print(dict_rNo2rName)

	### read mass from ff_file
	try:
		parse = ff_file.split(" ")
	except:
		nu.die("Error occurred when reading the force field file.. Check your " + str(ff_file))
	else:
		if not silent: print("Found " + str(len(parse)) + " Cerius2 Force Fields.")

	for i in parse:
		FF = dreiding.loadFF(i)
		temp_atominfo = dreiding.loadAtomTypes(FF)
		atominfo.update(temp_atominfo)


	### read trajectory file
	# how many steps to go?
	wc_trj_file = popen("grep TIMESTEP " + trj_file + " | wc -l ").read()
	n_timestep = int(wc_trj_file.split()[0]);
	print("The trajectory contains " + str(n_timestep) + " timesteps.")

	# Find header of the trajectory file
	while 1:
		templine = myTRJ.readline()
		line.append(templine.strip('\n').strip('ITEM: '))
		n_header += 1
		if "ITEM: ATOMS" in templine:
			break;

	# INITIAL trajectory information
	timestep = int(line[1])
	natoms = int(line[3])
	boxsize = [line[5].split(' ')[0], line[5].split(' ')[1], line[6].split(' ')[0], line[6].split(' ')[1], line[7].split(' ')[0], line[7].split(' ')[1]]
	boxsize = [float(i) for i in boxsize]
	keywords = line[8].strip('ATOMS ')

	# for every shot in the trajectory file update BGF and manipulate
	dumpatom = get_line(trj_file)
	processed_step = 0;

	t1 = t2 = 0; elapsed_time = 0;


	while 1:
		try:
			chunk = [next(dumpatom) for i in range(natoms+n_header)]
		except StopIteration:
			break;
	
		timestep = int(chunk[1])
		l_timestep.append(timestep)
		natoms = int(chunk[3])
		boxsize = [chunk[5].split(' ')[0], chunk[5].split(' ')[1], chunk[6].split(' ')[0], chunk[6].split(' ')[1], chunk[7].split(' ')[0], chunk[7].split(' ')[1]]; 
		boxsize = [float(i) for i in boxsize]; boxsize = [(boxsize[1] - boxsize[0]), (boxsize[3] - boxsize[2]), (boxsize[5] - boxsize[4])]
		keywords = chunk[8].split('ATOMS ')[1].strip('\n').split(' ')

		########
		if timestep < avg_timestep:
			continue;
		if timestep % 10000 != 0:
			continue;
		########

		### Show progress
		t1 = time.time();
		remaining_time = elapsed_time * (n_timestep - processed_step)
		sys.stdout.write('\r' + "Reading timestep.. " + str(timestep) + " (Elapsed time for the previous step: " + "{0:4.1f}".format(elapsed_time) + " seconds, Remaining time: " + "{0:4.1f}".format(remaining_time) + " seconds = " + "{0:4.1f} minutes".format(remaining_time/60) + ")")
		sys.stdout.flush()

		processed_step += 1;

		### update myBGF with trajectory information ###
		natom_bgf = len(myBGF.a)	# number of atoms in BGF file
	
		if not natom_bgf == natoms:
			nu.die("Number of atoms in trajectory file does not match with BGF file.")
	
		mode = ""
		if 'xs' in keywords or 'ys' in keywords or 'zs' in keywords:
			mode = 'scaled'
		elif 'x' in keywords or 'y' in keywords or 'z' in keywords:
			mode = 'normal'
		elif 'xu' in keywords or 'yu' in keywords or 'zu' in keywords:
			mode = 'unwrapped'
	
		# actual coordinate
		coordinfo = chunk[9:]
	
		# assume that coordinfo is similar to ['id', 'type', 'xs', 'ys', 'zs', 'ix', 'iy', 'iz']
		for atomline in coordinfo:
			atomcoord = atomline.split(' ')
			atom = myBGF.getAtom(int(atomcoord[0]))
	
			if mode == 'scaled':
				atom.x = float(atomcoord[2]) * boxsize[0]
				atom.y = float(atomcoord[3]) * boxsize[1]
				atom.z = float(atomcoord[4]) * boxsize[2]
	
			elif mode == 'unwrapped':
				atom.x = float(atomcoord[2])
				atom.y = float(atomcoord[3])
				atom.z = float(atomcoord[4])
	
			elif mode == 'normal':
				try:
					ix_index = keywords.index('ix')
					iy_index = keywords.index('iy')
					iz_index = keywords.index('iz')
				except ValueError:
					nu.warn("No image information no the trajectory file. Will be treated as unwrapped.")
					atom.x = float(atomcoord[2])
					atom.y = float(atomcoord[3])
					atom.z = float(atomcoord[4])
				else:
					atom.x = (int(atomcoord[ix_index]) * boxsize[0]) + float(atomcoord[2]) 
					atom.y = (int(atomcoord[iy_index]) * boxsize[1]) + float(atomcoord[3]) 
					atom.z = (int(atomcoord[iz_index]) * boxsize[2]) + float(atomcoord[4]) 
	
			try:
				for i in range(0, 3):
					myBGF.CRYSTX[i] = boxsize[i]
			except:
				pass;
				#nu.warn("Crystal information error: is this file not periodic?")
				
		# apply periodic condition
		myBGF = bgftools.periodicMoleculeSort(myBGF, 0)

		### DONE for updating trj on BGF file ###


		### find bin  ## CAUTION: only for NVT
		min = 100000; max = -100000
		for atom in myBGF.a:
			coord = [atom.x, atom.y, atom.z];
			if coord[axis] < min:
				min = coord[axis]
			if coord[axis] > max:
				max = coord[axis]

		if boxsize[axis] > max:
			max = boxsize[axis]
		if 0 < min:
			min = 0

		bins = np.arange(math.floor(min), math.ceil(max), interval)

		#print("Bins....:")
		#print(bins)


		### find CM of every molecule
		res_z = [];
		residues = [];
		for molecule in molecules:
			n_atom = len(molecule)	# number of atoms in the molecule
			residue_no = myBGF.getAtom(molecule[0]).rNo;	# residue number of the molecule
			m = bgftools.getMass(myBGF, molecule, ff_file);	# molecule mass
			cm = [0, 0, 0];		# center of mass
			coord = [0, 0, 0];	# atom x, y, z
			residues.append(residue_no)			# store residue numbers of molecules

			for aNo in molecule:
				atom = myBGF.getAtom(aNo)
				coord = [atom.x, atom.y, atom.z];

				# error if residue numbers are different
				if atom.rNo != residue_no:
					nu.die("Residue numbers in a same molecule are not consistent. Please check your BGF file.")

				# error if mass is unreadable
				ffType = string.strip(atom.ffType)
				try:
					aMass = atominfo[ffType]['MASS']
				except:
					nu.die("Cannot read the atom type " + str(atom.ffType) + " in the data file.")

				# calculate CM
				for index, i in enumerate(coord):
					cm[index] += i * aMass

			cm = [i / m for i in cm];	# center of mass

			res_z.append([residue_no, cm[2]])
		

		### binning according to z-axis
		bin_resno = [];
		for index, bin in enumerate(bins):
			temp_resno = [];
			for i in res_z:
				try:
					if bin < i[1] and bins[index+1] > i[1]:	
						temp_resno.append(i[0])
				except:
					continue
			bin_resno.append([bin, temp_resno])


		### refinement for residues
		bin_avg = [];
		for r in residue:
			temp_per_r = [];
			for b in bin_resno:
				temp_bin_res = [];
				for i in b[1]:
					if dict_rNo2rName[i] == r:
						temp_bin_res.append(i)
				temp_per_r.append([b[0], temp_bin_res])
			bin_avg.append([r, temp_per_r])


		### average
		#result = [];
		result_t = dict();
		for rdata in bin_avg:
			temp_per_r = [];
			for b in rdata[1]:
				temp_bin_res = [];
				avg_mag = 0;
				avg_angle_cos = 0;

				for i in b[1]:
					avg_mag += Dipole[timestep][rdata[0]]['ABSMU'][i]

					mu = Dipole[timestep][rdata[0]]['MU'][i]
					# angle btwn z-axis and mu
					angle_cos = mu[2] / math.sqrt(mu[0]**2 + mu[1]**2 + mu[2]**2)
					f_debug.write(str(timestep) + '\t' + str(mu) + '\t' + str(Dipole[timestep][rdata[0]]['ABSMU'][i]) + '\t' + str(angle_cos) + '\n')
					avg_angle_cos += angle_cos
				if not len(b[1]) == 0: avg_mag /= len(b[1])
				if not len(b[1]) == 0: avg_angle_cos /= len(b[1])
				temp_per_r.append([b[0], avg_angle_cos, avg_mag])	
			#result.append(temp_per_r)
			result_t[rdata[0]] = temp_per_r

		result[timestep] = result_t


		### end of loop: check elapsed time
		t2 = time.time()	# time mark
		elapsed_time = t2 - t1;


	### write pickle
	o = open(out_file + ".pickle", 'w')
	pickle.dump(result, o)

	### return

	print('')
	return 1

	### end of function


if __name__ == "__main__":
	bgf_file = ""; trj_file = ""; ff_file = ""; out_file = ""; pickle_file = ""; interval = 1; avg_timestep = 0;

	usage = """
	BGF_dipole.py -b bgf_file -t trj_file -f ff_file -o out_file -p pickle 
	"""

	options, args = getopt.getopt(sys.argv[1:], 'hb:t:f:o:p:a:', ['help', 'bgf=', 'trj=', 'ff=', 'out=', 'pickle=', 'avg='])

	if len(sys.argv) < 2:
		print(usage)
		sys.exit(0)

	print "Requested options: " + str(options)

	for option, value in options:
	        if option in ('-h', '--help'):
	                print(usage)
			sys.exit(0)
	        elif option in ('-b', '--bgf'):
	                bgf_file = value
	        elif option in ('-t', '--trj'):
	                trj_file = value
	        elif option in ('-f', '--ff'):
	                ff_file = value
	        elif option in ('-o', '--out'):
	                out_file = value
		elif option in ('-p', '--pickle'):
			pickle_file = value
		elif option in ('-a', '--avg'):
			avg_timestep = int(value)
	        elif option == NULL:
			print(usage)
			sys.exit(0)

	if out_file == "":
		out_file = bgf_file[:-4] + ".dipoleDistr";
	
	# main call
	densityProfile(bgf_file, trj_file, ff_file, pickle_file, out_file, avg_timestep)
