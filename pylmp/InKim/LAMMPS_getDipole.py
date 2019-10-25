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
import cPickle as pkl

import dump

option = ""; args = ""; bgf_file = ""; trj_file = ""; out_file = ""; 
usage = """
BGF_dipole.py -b bgf_file -t trj_file -f ff_file
"""
version = "130206"

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


def dipole(bgf_file, trj_file, ff_file, out_file, avg_timestep, silent=False):

	### init
	timestep = 0; l_timestep = []; line = []; n_header = 0;
	t1 = 0; t2 = 0; # clock
	vector = [0, 0, 1];	# the axis of interest. in the acetone-water case, z direction.
	atominfo = dict();	# atom data extracted from ff_file
	result = dict();
	axis = 2;		# for orientational distribution of dipole moments. 1: x-axis, 2: y-axis, 3: z-axis


	### unit conversion for debye
	elementary_q = 1.602176487e-19 # elementary charge in C
	debye_conv = 3.33564e-30       # 1 Debye in C*m
	k = elementary_q*1e-10/debye_conv


	### open files
	myBGF = bgf.BgfFile(bgf_file)
	myTRJ = open(trj_file)
	myTRJ.seek(0)
	myOUT = open(out_file + ".dat", 'w')
	#myOUT2 = open(out_file + ".angle.dat", 'w')
	myRESULT = open(out_file + ".pickle", 'w')
	myRESULT2 = open(out_file + ".histo.pickle", 'w')


	### read residues from bgf_file
	residueNames = set();	# kind of residues in BGF file
	residueNumbers = set();
	dict_residue = dict();	# stores residue numbers per each residue. RES1: [1, 2, ..], RES2: [6, 7, ..]
	for i in myBGF.a:
		rname = string.strip(i.rName)
		residueNames.add(rname)
		
		rno = i.rNo
		residueNumbers.add(rno)

	temp = "";
	for i in residueNames:
		temp += i + " ";
		dict_residue[i] = [];

	if not silent: print("Found " + str(len(residueNames)) + " residues in BGF file: " + str(temp))


	### bookkeep residue numbers for each residue
	molecules = bgftools.getMoleculeList(myBGF)
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
	wc_trj_file = popen("grep -A 1 TIMESTEP " + trj_file).read()
	wc_trj_file = wc_trj_file.split()
	l_timestep = [];

	for i in wc_trj_file:
		if "ITEM" in i or "TIME" in i or "--" in i:
			pass;
		else:
			l_timestep.append(i)

	l_timestep = [ int(i) for i in l_timestep ]
	l_timestep.sort()

	n_timestep = len(l_timestep)
	if not silent: print("The trajectory contains " + str(n_timestep) + " timesteps.")

	l_requested_timesteps = l_timestep[-avg_timestep:]	# requested timesteps
	if not silent: print("Only requested the last " + str(avg_timestep) + " timesteps. ")

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

		if not timestep in l_requested_timesteps:
			continue;

		#l_timestep.append(timestep)
		natoms = int(chunk[3])
		boxsize = [chunk[5].split(' ')[0], chunk[5].split(' ')[1], chunk[6].split(' ')[0], chunk[6].split(' ')[1], chunk[7].split(' ')[0], chunk[7].split(' ')[1]]; 
		boxsize = [float(i) for i in boxsize]; boxsize = [(boxsize[1] - boxsize[0]), (boxsize[3] - boxsize[2]), (boxsize[5] - boxsize[4])]
		keywords = chunk[8].split('ATOMS ')[1].strip('\n').split(' ')

		### Show progress
		t1 = time.time();
		remaining_time = elapsed_time * (len(l_requested_timesteps) - processed_step)
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


		### Now Dipole Moment Calculation
		# list of all molecules in BGF
		dipole_moments = dict();
		abs_dipole_moments = dict();
		angles = dict();
		temp_dict = dict();
		output_dipole = str(timestep) + "\t";
		output_abs_dipole = str(timestep) + "\t";
		output_angle = str(timestep) + "\t";


		### find bin for orientational distribution of dipole moment
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

		bins = np.arange(math.floor(min), math.ceil(max), 1.0)
		binned_rNo = dict();
		for r in residueNames:
			binned_rNo[r] = [ [] for b in bins];	# space for binned residue numbers (= molecules)
		binned_mu = copy.deepcopy(binned_rNo)
		

		### For every molecule
		for molecule in molecules:
			n_atom = len(molecule)	# number of atoms in the molecule
			residue_no = myBGF.getAtom(molecule[0]).rNo;	# residue number of the molecule
			residue_name = myBGF.getAtom(molecule[0]).rName;	# residue name of the molecule
			m = bgftools.getMass(myBGF, molecule, ff_file);	# molecule mass
			cm = [0, 0, 0];		# center of mass
			mu = [0, 0, 0];		# dipole moment
			coord = [0, 0, 0];	# atom x, y, z

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
			#print(str(residue_no) + '\t' + str(cm))

			# REMARK: it seems to be okay without this process.
			"""
			# translate the molecule to CM
			for aNo in molecule:
				atom = myBGF.getAtom(aNo)
				atom.x -= cm[0]
				atom.y -= cm[1]
				atom.z -= cm[2]
			"""

			# calculate dipole moment sum(qi * di)
			mu_i = [];
			for aNo in molecule:
				atom = myBGF.getAtom(aNo)
				coord = [atom.x, atom.y, atom.z];
				mu_i.append([i * atom.charge * k for i in coord]);

			# sum up
			mu = [0, 0, 0];
			for i in mu_i:
				mu[0] += i[0]
				mu[1] += i[1]
				mu[2] += i[2]


			# for acetone case
			#len_mu = math.sqrt(dot(mu, mu))
			#len_vector = math.sqrt(dot(vector, vector))
			len_mu = math.sqrt(dot(mu, mu))
			len_vector = 1.0;

			#angle = dot(mu, vector) / ( len_mu * len_vector ) 
			angle = mu[2] / (len_mu * len_vector)	# for ACETONE case

			# results
			dipole_moments[residue_no] = mu;
			abs_dipole_moments[residue_no] = len_mu;
			angles[residue_no] = angle;


			### Orientational distribution of dipole moments
			# binning according to molecule's CM and store to binned_rNo
			for index, bin in enumerate(bins):
				try:
					if bin < cm[axis] and bins[index+1] > cm[axis]:	
						binned_rNo[residue_name][index].append(residue_no);
				except:
					continue;
			# binning is perfect so far!


		residueNumbers = list(residueNumbers)
		residueNames = list(residueNames)
		residueNumbers.sort()

		temp_dict = dict();
		temp_dict2 = dict();
		temp_dict2['MU'] = dipole_moments
		temp_dict2['ABSMU'] = abs_dipole_moments
		temp_dict2['ANGLES'] = angles
		temp_dict['TOTAL'] = temp_dict2

		### paperworks per residue
		for resname in residueNames:
			temp_dict2 = dict();
			res_dipole_moments = dict();
			res_abs_dipole_moments = dict();
			res_angles = dict();

			# load residue numbers per residue type
			resnumbers = dict_residue[resname]
			for rno in resnumbers:
				res_dipole_moments[rno] = dipole_moments[rno]
				res_abs_dipole_moments[rno] = abs_dipole_moments[rno]
				res_angles[rno] = angles[rno]

			temp_dict2['MU'] = res_dipole_moments
			temp_dict2['ABSMU'] = res_abs_dipole_moments
			temp_dict2['ANGLES'] = res_angles
			temp_dict2['BIN'] = bins

			# for analyze
			#temp_dict2['RNO_DISTR'] = binned_rNo[resname]
			temp_dict2['DISTR_ANGLES'] = binned_mu[resname]
			temp_dict2['DISTR_ABSMU'] = binned_mu[resname]

			# write dipole moments
			for index, item in enumerate(binned_rNo[resname]):
				if len(item) == 0:
					continue;
				for i in item:
					temp_dict2['DISTR_ANGLES'][index].append(temp_dict2['ANGLES'][i])
					#temp_dict2['DISTR_ABSMU'][index].append(temp_dict2['ABSMU'][i])

			temp_dict[resname] = temp_dict2

		# result: timestep - residue - MU, ANGLES, NATOMS
		result[timestep] = temp_dict


		t2 = time.time()	# time mark
		elapsed_time = t2 - t1;


	### Averaging orientational distribution of dipole moments over timesteps
	#del temp_dict1
	del temp_dict2
	temp_dict1 = dict();	# angles with keys: bin
	temp_dict2 = dict();	# |mu| with keys: bin

	# append
	for r in residueNames:
		temp_dict1[r] = dict()
		temp_dict2[r] = dict()

		for t in l_requested_timesteps:
			for index, b in enumerate(result[t][r]['BIN']):
				if not temp_dict1[r].has_key(b):
					temp_dict1[r][b] = [];
				if not temp_dict2[r].has_key(b):
					temp_dict2[r][b] = [];

				# append averaged values
				avg, std = meanstdv(result[t][r]['DISTR_ANGLES'][index])
				temp_dict1[r][b].append(avg)

				for i in result[t][r]['DISTR_ANGLES'][index]:
					temp_dict2[r][b].append(i)

				#for i in result[t][r]['DISTR_ABSMU'][index]:
				#	temp_dict2[b].append(i)
	
	# average
	avg_angles = dict();
	avg_absmu = dict();
	avg_angles_stdev = dict();
	avg_absmu_stdev = dict();

	for r in residueNames:
		avg_angles[r] = dict();
		avg_absmu[r] = dict();
		avg_angles_stdev[r] = dict();
		avg_absmu_stdev[r] = dict();
		temp = temp_dict1[r].keys()
		temp.sort()
		for i in temp:
			avg_angles[r][i], avg_angles_stdev[r][i] = meanstdv(temp_dict1[r][i])
			avg_absmu[r][i], avg_absmu_stdev[r][i] = meanstdv(temp_dict2[r][i])

	# print for average
	print("")
	print("Averaged " + str(avg_timestep) + " timesteps out of " + str(n_timestep))
	for r in residueNames:
		print(r)
		for i in temp:
			#print(str(i) + '\t' + str(temp_dict1[i]))
			print(str(i) + '\t' + str(avg_angles[r][i]) + '\t' + str(avg_angles_stdev[r][i]))
		print("")


	# save for debug
	output = ""
	for r in residueNames:
		for i in temp:
			output += str(i) + '\t' + str(temp_dict1[r][i]) + '\n'
	
	myOUT.write(output)
	myOUT.close()


	# save the pickle object
	pkl.dump(result, myRESULT)
	myRESULT.close()

	pkl.dump(temp_dict2, myRESULT2)
	myRESULT2.close()


	print('')
	return 1

	### end of function

def meanstdv(x):
	from math import sqrt
	n, mean, std = len(x), 0.0, 0.0

	if n == 0:
		return 0.0, 0.0

	for a in x:
		mean = mean + a
	mean = mean / float(n)

	for a in x:
		std = std + (a - mean)**2

	if n != 1:
		std = sqrt(std / float(n-1))
	else:
		std = 0.0

	return mean, std


def dot(a, b):
	"""
	Calculate dot product of a and b
	"""
	try:
		na = len(a); nb = len(b)
	except:
		nu.die("List required for dot product calculation.")

	na = len(a); nb = len(b)

	temp = [a[i] * b[i] for i in range(na)]
	temp2 = 0;

	for i in temp:
		temp2 += i

	return temp2;


def rad2deg(r):
	return 180*(r/math.pi)


if __name__ == "__main__":
	bgf_file = ""; trj_file = ""; ff_file = ""; out_file = ""; avg_timestep = 0;

	options, args = getopt.getopt(sys.argv[1:], 'hb:t:f:o:n:', ['help', 'bgf=', 'trj=', 'ff=', 'out=', 'average='])

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
		elif option in ('-n', '--avg'):
			avg_timestep = int(value)
	        elif option == NULL:
			print(usage)
			sys.exit(0)

	if out_file == "":
		out_file = bgf_file[:-4] + ".dipole";
	if avg_timestep == 0:
		print("The result will be averaged over all timesteps")
	
	# main call
	dipole(bgf_file, trj_file, ff_file, out_file, avg_timestep)
