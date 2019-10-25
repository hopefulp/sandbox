#!/home/noische/program/python27/bin/python

import sys, re, string, getopt, optparse, math, time, pprint 
from os import popen
import bgf
import bgftools
import numpy
import operator
import copy
import nutils as nu
import itertools
import time

import dump

option = ""; args = ""; bgf_file = ""; trj_file = ""; out_file = ""; 
usage = """
countWaterCNT.py -b bgf_file -t trj_file -f ff_file(s) -n #step -H do_hbond
"""
version = "130123"

#-----------------
# calculate the distance of specified atoms in every lammps trajectory file
# this is only used during pei analysis by in kim 
# last update: 2012/06/01
#_________________
def get_line(trj_file):
	"""
	Open a lazy file stream. This is considered to increase a file reading speed.
	"""
	with open(trj_file) as file:
		for i in file:
			yield i

def sortkey(list):
	return list[0];

def countWaterCNT(bgf_file, trj_file, ff_file, n_step, do_hbonds, silent=False):

	# const
	PI = math.pi
	vdw_r_C = 1.7

	# init
	timestep = 0; l_timestep = []; line = []; n_header = 0;
	t1 = 0; t2 = 0; # clock

	myBGF = bgf.BgfFile(bgf_file)
	myTRJ = open(trj_file)
	myTRJ.seek(0)
	myDAT = open(bgf_file[:-4] + ".calcDensity.count.dat", "w")
	myDAT2 = open(bgf_file[:-4] + ".calcDensity.hbond.dat", "w")
	myDAT3 = open(bgf_file[:-4] + ".calcDensity.configuration.dat", "w")
	myDAT.write(str(sys.argv) + "\n")
	myDAT.write("Time\tmin_X\tmax_X\tHeight\tRadius\tNumber\tVolume\tEff_volume\tExact_mass\tMass\tExact_den\tMtOH_den\tExact_eff_den\tMtOH_eff_den\tn_hbond\n")

	# how many steps to go?
	wc_trj_file = popen("grep TIMESTEP " + trj_file + " | wc -l ").read()
	n_timestep = int(wc_trj_file.split()[0]);
	if n_step == 0:
		n_step = n_timestep;

	print("The trajectory contains " + str(n_timestep) + " timesteps.")
	print("The script will proceed for the last " + str(n_step) + " timesteps.")

	# extract aNos of CNT in the BGF file
	aNo_CNT = []
	aNo_MtOH_C = []
	aNo_MtOH_all = []

	for atom in myBGF.a:
		# Carbons in CNT
		if "CNT" in atom.rName:
			aNo_CNT.append(atom.aNo)

		# Carbons in MeOH
		if "MET" in atom.rName and "C" in atom.aName:
			aNo_MtOH_C.append(atom.aNo)

	N_CNT = len(aNo_CNT)	# the number of CNT atoms

	# check if there exists MtOH properly
	if len(aNo_MtOH_C) == 0:
		nu.die("No MtOH molecules in the BGF file.")
	if len(aNo_CNT) == 0:
		nu.die("No CNT molecules in the BGF file.")

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
		### Show progress
		t1 = time.time();
		remaining_time = elapsed_time * (n_step - processed_step)
		sys.stdout.write('\r' + "Reading timestep.. " + str(timestep) + " (Remaining time: " + "{0:4.1f}".format(remaining_time) + " seconds, " + "{0:4.1f} minutes".format(remaining_time/60) + ")")
		sys.stdout.flush()

		processed_step += 1;
		if processed_step == n_step:
			break;

		aNo_MtOH_C_atoms = []
		aNo_atoms_in_CNT = []


		### Read
		try:
			chunk = [next(dumpatom) for i in range(natoms+n_header)]
		except StopIteration:
			break;

		timestep = int(chunk[1])
		natoms = int(chunk[3])
		boxsize = [chunk[5].split(' ')[0], chunk[5].split(' ')[1], chunk[6].split(' ')[0], chunk[6].split(' ')[1], chunk[7].split(' ')[0], chunk[7].split(' ')[1]]; 
		boxsize = [float(i) for i in boxsize]; boxsize = [(boxsize[1] - boxsize[0]), (boxsize[3] - boxsize[2]), (boxsize[5] - boxsize[4])]
		keywords = chunk[8].split('ATOMS ')[1].strip('\n').split(' ')


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
		myBGF = bgftools.periodicMoleculeSort(myBGF, 0, True)

		### myBGF update complete! ###

		aNo_MtOH_C_in_CNT = copy.deepcopy(aNo_MtOH_C)

		### Realign whole system to x axis
		sys.stdout.write(' Realigning.. ')
		sys.stdout.flush()
		# initialize for the moment of inertia and the center of mass calculation
		U = 0; Ut = 0; Uv = 0;
		Ixx = 0; Ixy = 0; Ixz = 0;
		Iyx = 0; Iyy = 0; Iyz = 0;
		Izx = 0; Izy = 0; Izz = 0;
		Mx = 0; My = 0; Mz = 0;

		# transpose for "all atoms in BGF": move COM of CNT as origin
		# com of CNT
		Mx = 0; My = 0; Mz = 0;
		for atom in myBGF.a:
			if "CNT" in atom.rName:
				Mx += atom.x / N_CNT
				My += atom.y / N_CNT
				Mz += atom.z / N_CNT

		# move
		for atom in myBGF.a:
			atom.x -= Mx
			atom.y -= My
			atom.z -= Mz

		# calculate the moment of inertia (MI) and the center of mass (COM) of CNT from "myBGF"
		for aNo in aNo_CNT:
			atom = myBGF.getAtom(aNo)
			Ixx += (atom.y**2 + atom.z**2) / N_CNT
			Iyy += (atom.x**2 + atom.z**2) / N_CNT
			Izz += (atom.x**2 + atom.y**2) / N_CNT
			Ixy -= (atom.x * atom.y) / N_CNT
			Ixz -= (atom.x * atom.z) / N_CNT
			Iyz -= (atom.y * atom.z) / N_CNT
			
		# the moment of inertia tensor
		I = numpy.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])

		# eigenvalue & eigenvector calculation
		eigval, eigvec = numpy.linalg.eig(I)	# eigval[0] is the minimum among the values.
	
		# rearrange the U vector
		U = numpy.matrix(eigvec)
		Ut = U.T

		# "myBGF" rotation
		for atom in myBGF.a:
			v = numpy.matrix([atom.x, atom.y, atom.z]).T
			Uv = Ut * v
			atom.x = float(Uv[0])
			atom.y = float(Uv[1])
			atom.z = float(Uv[2])

		# for CNT atoms, calculate some properties
		min_x_CNT = 0.0; max_x_CNT = 0.0; radius_CNT = 0.0; height_CNT = 0.0; aNo_MtOH_C_not_in_CNT = []; temp = [];
		min_y_CNT = 0.0; max_y_CNT = 0.0;
		min_z_CNT = 0.0; max_z_CNT = 0.0;
		CNT_orientation = "";

		# check the orientation of CNT
		for aNo in aNo_CNT:
			atom = myBGF.getAtom(aNo)
			# minimum and maximum x coord of CNT: this will be the height of CNT
			if atom.x < min_x_CNT:
				min_x_CNT = atom.x
			if atom.x > max_x_CNT:
				max_x_CNT = atom.x
			if atom.y < min_y_CNT:
				min_y_CNT = atom.y
			if atom.y > max_y_CNT:
				max_y_CNT = atom.y
			if atom.z < min_z_CNT:
				min_z_CNT = atom.z
			if atom.z > max_z_CNT:
				max_z_CNT = atom.z
	
		x_diff = max_x_CNT - min_x_CNT
		y_diff = max_y_CNT - min_y_CNT
		z_diff = max_z_CNT - min_z_CNT

		if x_diff > y_diff and x_diff > z_diff:
			# CNT is aligned along x axis
			height_CNT = x_diff
			CNT_orientation = "x"
		elif y_diff > x_diff and y_diff > z_diff:
			# CNT is aligned along y axis
			height_CNT = y_diff
			CNT_orientation = "y"
		elif z_diff > x_diff and z_diff > y_diff:
			# CNT is aligned along z axis
			height_CNT = z_diff
			CNT_orientation = "z"

		if CNT_orientation == "x":
			for aNo in aNo_CNT:
				atom = myBGF.getAtom(aNo)
				radius_CNT += math.sqrt(atom.y**2 + atom.z**2)	# average radius of CNT
		elif CNT_orientation == "y":
			for aNo in aNo_CNT:
				atom = myBGF.getAtom(aNo)
				radius_CNT += math.sqrt(atom.x**2 + atom.z**2)	# average radius of CNT
		elif CNT_orientation == "z":
			for aNo in aNo_CNT:
				atom = myBGF.getAtom(aNo)
				radius_CNT += math.sqrt(atom.x**2 + atom.y**2)	# average radius of CNT

		radius_CNT = radius_CNT / N_CNT

		### Rotate y or z to x axis
		Tyx = numpy.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
		Tzx = numpy.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
	
		if CNT_orientation == "y":
			for atom in myBGF.a:
				u = numpy.matrix([atom.x, atom.y, atom.z]).T
				Tu = Tyx*u
				atom.x = float(Tu[0])
				atom.y = float(Tu[1])
				atom.z = float(Tu[2])
				min_x_CNT = min_y_CNT;
				max_x_CNT = max_y_CNT;
	
		elif CNT_orientation == "z":
			for atom in myBGF.a:
				u = numpy.matrix([atom.x, atom.y, atom.z]).T
				Tu = Tzx*u
				atom.x = float(Tu[0])
				atom.y = float(Tu[1])
				atom.z = float(Tu[2])
				min_x_CNT = min_z_CNT;
				max_x_CNT = max_z_CNT;
	
		### get MtOHs in CNT
		sys.stdout.write('Density.. ')
		sys.stdout.flush()
		# 1. aNo_atoms_in_CNT: exact atoms within CNT
		for atom in myBGF.a:
			ano = atom.aNo
			dist = math.sqrt(atom.y**2 + atom.z**2)
			if not "CNT" in atom.rName and atom.x > min_x_CNT and atom.x < max_x_CNT and dist < radius_CNT:
				# inside
				aNo_atoms_in_CNT.append(ano)
			else:
				# outside
				pass;

		# 2. aNo_MtOH_C_atoms: molecules which C atom is within CNT
		for aNo in aNo_MtOH_C:
			atom = myBGF.getAtom(aNo)
			dist = math.sqrt(atom.y**2 + atom.z**2)
			if atom.x > min_x_CNT and atom.x < max_x_CNT and dist < radius_CNT:
				# inside of the CNT
				pass;
			else:
				"""
				# delete atoms
				delete_list = [];
				dummy = bgftools.getmolecule(myBGF, myBGF.getAtom(aNo), delete_list)
				for i in delete_list:
					aNo_MtOH_C_not_in_CNT.append(myBGF.a2i[i])
				"""
				# outside of the CNT
				aNo_MtOH_C_in_CNT.remove(aNo)

		for aNo in aNo_MtOH_C_in_CNT:
			anos = [];
			dummy = bgftools.getmolecule(myBGF, myBGF.getAtom(aNo), anos)
			for i in anos:
				aNo_MtOH_C_atoms.append(i)

		# calculate volume of CNT
		# radius: average radius - can't be applied in bent CNTs
		# height: maximum height of CNT
		# volume: radius**2 + height
		volume = (radius_CNT**2) * height_CNT * PI
		eff_volume = (radius_CNT - vdw_r_C)**2 * height_CNT * PI

		# calculate mass in CNT
		MtOH_exact_mass = bgftools.getMass(myBGF, aNo_atoms_in_CNT, ff_file)
		MtOH_MtOH_C_mass = bgftools.getMass(myBGF, aNo_MtOH_C_atoms, ff_file)

		# calculate density: mass / 6.022 / vol * 10
		vol_exact_density = MtOH_exact_mass / 6.022 / volume * 10;
		vol_MtOH_density = MtOH_MtOH_C_mass / 6.022 / volume * 10;

		# calculate effective density ( radius - r_vdw )
		eff_vol_exact_density = MtOH_exact_mass / 6.022 / eff_volume * 10;
		eff_vol_MtOH_density = MtOH_MtOH_C_mass / 6.022 / eff_volume * 10;

		### Check MET status
		output = "";
		MeOH_align = [];

		# find O in MET: MeOH_align = [C.x, O.x]
		for ano in aNo_MtOH_C_in_CNT:
			atom = myBGF.getAtom(ano)
			for i in atom.CONECT:
				atom2 = myBGF.getAtom(i)
				if "O" in atom2.ffType:
					MeOH_align.append([atom.x, atom2.x])
	
		MeOH_align.sort(key=sortkey)

		# print MET align status
		for i in MeOH_align:
			if i[0] > i[1]:
				# C-O: 1
				output += "1 "
			else:
				# O-C: 0
				output += "0 "

		# n_hbond
		n_hbond = 0;
		"""
		# how many 10 or 01 (hydrogen bonds) in output?
		count += output.count("1 0")
		count += output.count("0 1")
		"""

		hbond_C_pairs = "";
		hbond_OH_pairs = [];
		if do_hbonds:
			### find hydrogen bonds
			### criteria: r(O..O) <= 3.5, r(O..HO) <= 2.6, A(HO-O..O) <= 30' (ref: Haughney et al., JPC 1987)
			###                       acceptor donor       donor   acceptor
			# for all pairs of MtOH molecules
			sys.stdout.write('Hbonds.. ')
			sys.stdout.flush()

			C_pairs = itertools.permutations(aNo_MtOH_C_in_CNT, 2)
			# count how many iterations on the permutations:
			#c_perm = 0;
			#for i in C_pairs:
			#	c_perm += 1;

			#print(str(len(aNo_MtOH_C_in_CNT)) + " " + str(c_perm))

			for p in C_pairs:
				(ano1, ano2) = p

				# Donor part
				# get C1
				C1 = myBGF.getAtom(ano1)
				# get O1
				for i in C1.CONECT:
					temp_atom = myBGF.getAtom(i)
					if "O" in temp_atom.ffType[0]:
						O1 = temp_atom
				# get H1
				for i in O1.CONECT:
					temp_atom = myBGF.getAtom(i)
					if "H" in temp_atom.ffType[0]:
						H1 = temp_atom

				# Acceptor part
				# get C2
				C2 = myBGF.getAtom(ano2)
				# get O2
				for i in C2.CONECT:
					temp_atom = myBGF.getAtom(i)
					if "O" in temp_atom.ffType[0]:
						O2 = temp_atom
				# get H2
				for i in O2.CONECT:
					temp_atom = myBGF.getAtom(i)
					if "H" in temp_atom.ffType[0]:
						H2 = temp_atom

				"""
				# H1 and O2 should be inside the CNT
				dist = math.sqrt(H1.y**2 + H1.z**2)
				if H1.x > min_x_CNT and H1.x < max_x_CNT and dist < radius_CNT:
					continue;

				dist = math.sqrt(O2.y**2 + O2.z**2)
				if O2.x > min_x_CNT and O2.x < max_x_CNT and dist < radius_CNT:
					continue;
				"""

				## determine hbond
				# calculate r(O1..O2)
				crit_1 = bgf.distance(O1, O2);
				# calculate r(O1..H2)
				crit_2 = bgf.distance(O2, H1);
				# calculate A(O1-H2..O2)
				crit_3 = bgf.angle(H1, O1, O2);

				#if crit_1 <= 3.5 and crit_2 <= 2.6 and crit_3 <= 30:
				if crit_1 <= 3.5 and crit_3 <= 30:
					# Hbond!!
					n_hbond += 1;
					OH_pair = [O2.aNo, H1.aNo];
					hbond_OH_pairs.append(OH_pair)
					#hbond_C_pairs += str(ano1) + "-" + str(ano2) + "\t"
				else:
					continue;

		hbond_pairs = hbond_OH_pairs
		#temp = nu.removeRepeat(hbond_OH_pairs)
		#hbond_pairs = nu.removeReverse(temp)

		sys.stdout.write('Done                                                        ')
		sys.stdout.flush()

		myDAT.write(str(timestep) + "\t")
		myDAT.write(str("{0:<8.3f}".format(min_x_CNT)) + "\t")
		myDAT.write(str("{0:<8.3f}".format(max_x_CNT)) + "\t")
		myDAT.write(str("{0:<8.3f}".format(height_CNT)) + "\t")
		myDAT.write(str("{0:<8.3f}".format(radius_CNT)) + "\t")
		myDAT.write(str(len(aNo_MtOH_C_in_CNT)) + "\t")
		myDAT.write(str("{0:<8.3f}".format(volume)) + "\t")
		myDAT.write(str("{0:<8.3f}".format(eff_volume)) + "\t")
		myDAT.write(str("{0:<8.3f}".format(MtOH_exact_mass)) + "\t")
		myDAT.write(str("{0:<8.3f}".format(MtOH_MtOH_C_mass)) + "\t")
		myDAT.write(str("{0:<8.3f}".format(vol_exact_density)) + "\t")
		myDAT.write(str("{0:<8.3f}".format(vol_MtOH_density)) + "\t")
		myDAT.write(str("{0:<8.3f}".format(eff_vol_exact_density)) + "\t")
		myDAT.write(str("{0:<8.3f}".format(eff_vol_MtOH_density) + "\t"))
		myDAT.write(str(len(hbond_pairs)) + "\n")
		myDAT2.write(str(hbond_pairs) + "\n")
		myDAT3.write(str(output) + "\n")
		#myDAT.write(str(aNo_MtOH_C_in_CNT) + "\n")

		# write output
		#myBGF.saveBGF(bgf_file[:-4] + "." + str(timestep) + ".bgf")
		#myBGF2.saveBGF(bgf_file[:-4] + "." + str(timestep) + ".bgf")

		t2 = time.time()	# time mark
		elapsed_time = t2 - t1;

	print('')

	myDAT.close()
	myDAT2.close()

	return 1

	### end of function


if __name__ == "__main__":
	bgf_file = ""; trj_file = ""; ff_file = ""; n_step = 0; do_hbonds = False;

	options, args = getopt.getopt(sys.argv[1:], 'hb:t:f:n:H', ['help', 'bgf=', 'trj=', 'forcefield=', 'step=', 'hbonds'])

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
	        elif option in ('-f', '--forcefield'):
	                ff_file = str(value).strip()
		elif option in ('-n', '--step'):
			n_step = int(value)
		elif option in ('-H', '--hbonds'):
			do_hbonds = True
	        elif option == NULL:
			print(usage)
			sys.exit(0)
	
	# more options: resname CNT, resname solvent, save BGF

	# main call
	countWaterCNT(bgf_file, trj_file, ff_file, n_step, do_hbonds)
