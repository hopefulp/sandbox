#!/home/noische/python -u

import sys, re, string, getopt, optparse, math, time
from os import popen
import time
import os

import numpy as np

import bgf
import bgftools
import nutils as nu
import lammpstools as lt
import dreiding

option = ""; args = ""; bgf_file = ""; trj_file = ""; out_file = ""; 
usage = """
Count the number of molecules in CNT with LAMMPS trajectory file.

countWaterCNT.py -b bgf_file -t trj_file -n #step -C
	-b	BGF file with nanotube residue xNT and water residue WAT
	-t	LAMMPS trajectory file
	-n	(optional) Number of steps to calculate
	-C	Considers water movement within PBC. 
		** This takes a lot of time if the position of NT is not fixed. **
"""
version = "150102"

#-----------------
# Get averaged velocity of water molecules confined in CNT
# 
# 
#_________________
def get_line(trj_file):
	"""
	Open a lazy file stream. This is considered to increase a file reading speed.
	"""
	with open(trj_file) as file:
		for i in file:
			yield i


def countWater(bgf_file, trj_file, n_step, watercopy, ff_file, silent=False):

	### const
	PI = math.pi
	vdw_r_C = 1.7

	### init
	timestep = 0; l_timestep = []; line = []; n_header = 0;
	t1 = 0; t2 = 0; # clock

	myBGF = bgf.BgfFile(bgf_file)
	myTRJ = open(trj_file)
	myTRJ.seek(0)

	global out_file
	if out_file == "": out_file = "countWater5.profile"
	print("The result will be recorded to the file " + out_file + " ...")
	ftemp = open(out_file, 'w')
	ftemp.write(str(sys.argv) + "\n")
	ftemp.write("t" + '\t' + "n_O" + '\t' + "r_CNT" + '\t' + "std_r" + '\t' + "min_x" + '\t' + "max_x" + '\t' + "min_y" + '\t' + "max_y" + '\t' + "min_z" + '\t' + "max_z" + '\t' + "l_NT" + '\t' + "replNum" + '\t' + "copyNum" + '\t' + "n_WAT" + '\t' + "remark" + '\n')

	curr_dir = os.path.abspath(".")
	temp_dir = curr_dir + "/countWAT5/"
	print(temp_dir)
	if not os.path.isdir(temp_dir): os.makedirs(temp_dir)


	### how many steps to go?
	n_timestep = len(lt.getTrjInfo(trj_file))
	if n_step == 0:
		n_step = n_timestep;

	print(" ..The trajectory contains " + str(n_timestep) + " timesteps.")
	print("The script will proceed for the first " + str(n_step) + " timesteps.")


        ### read mass from ff_file
	atominfo = dict();
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


	### extract aNos of CNT in the BGF file
	aNo_CNT = []; aNo_WAT_O = []; aNo_WAT_all = []; mass_CNT = [];

	for atom in myBGF.a:
		# Carbons in CNT or atoms in BNNT
		if "NT" in atom.rName:
			aNo_CNT.append(atom.aNo)
			ffType = string.strip(atom.ffType)
			try:
				aMass = atominfo[ffType]['MASS']
			except:
				nu.die("Cannot read the atom type " + str(atom.ffType) + " in the data file.")
			mass_CNT.append(aMass)

		# Oxygen in water
		if "WAT" in atom.rName and "O" in atom.aName:
			aNo_WAT_O.append(atom.aNo)

	N_CNT = len(aNo_CNT)	# the number of CNT atoms


	### check if there exists water properly
	if len(aNo_WAT_O) == 0:
		nu.die("No water molecules in the BGF file.")
	if len(aNo_CNT) == 0:
		nu.die("No CNT molecules in the BGF file.")


	### Find header of the trajectory file
	while 1:
		templine = myTRJ.readline()
		line.append(templine.strip('\n').strip('ITEM: '))
		n_header += 1
		if "ITEM: ATOMS" in templine:
			break;


	### INITIAL trajectory information
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

		if processed_step == n_step:
			break;

		### Read
		myBGF = bgf.BgfFile(bgf_file)
		try:
			chunk = [next(dumpatom) for i in range(natoms+n_header)]
		except StopIteration:
			break;

		timestep = int(chunk[1])
		natoms = int(chunk[3])
		boxsize = [chunk[5].split(' ')[0], chunk[5].split(' ')[1], chunk[6].split(' ')[0], chunk[6].split(' ')[1], chunk[7].split(' ')[0], chunk[7].split(' ')[1]]; 
		boxsize = [float(i) for i in boxsize]; boxsize = [(boxsize[1] - boxsize[0]), (boxsize[3] - boxsize[2]), (boxsize[5] - boxsize[4])]
		keywords = chunk[8].split('ATOMS ')[1].strip('\n').split(' ')

		mode = 'unwrapped'	# assume that "dump            1 all custom 100 ${sname}${rtemp}K.nvt.lammps id type xu yu zu vx vy vz" in lammps input
	
		# actual coordinate
		coordinfo = chunk[9:]
	
		# load atom coordinates from chunk
		for atomline in coordinfo:
			atomcoord = atomline.split(' ')
			atom = myBGF.getAtom(int(atomcoord[0]))
	
			atom.x = float(atomcoord[2])
			atom.y = float(atomcoord[3])
			atom.z = float(atomcoord[4])

			
		### convert atom coordinates to lists
		CNT = []; WATER = []; 
		for atom in myBGF.a:
			if "NT" in atom.rName:
				CNT.append([atom.x, atom.y, atom.z])
			if "WAT" in atom.rName and "O" in atom.ffType:
				WATER.append([atom.x, atom.y, atom.z])

		CNT = np.array(CNT); n_CNT = len(CNT)	# CNT coordinates
		WATER = np.array(WATER)	# Water coordinates
		WATERONLY = np.copy(WATER)
		boxsize = np.array(boxsize)


		### initialize for the moment of inertia and the center of mass calculation
		U = 0; Ut = 0; Uv = 0;
		Ixx = 0; Ixy = 0; Ixz = 0;
		Iyx = 0; Iyy = 0; Iyz = 0;
		Izx = 0; Izy = 0; Izz = 0;
		Mx = 0; My = 0; Mz = 0;


		### transpose "all atoms in BGF": move COM of CNT as origin
		Mx, My, Mz = np.average(CNT, axis=0, weights=mass_CNT)	# CoM of CNT
		COM = np.array([[Mx, My, Mz]])

		CNT = CNT - COM + boxsize / 2.0	# move
		WATER = WATER - COM + boxsize / 2.0


		### apply PBC 
		WATER = np.mod(WATER, boxsize)


		### save coordinates to BGF before rotation
		myBGF2 = bgf.BgfFile()
		for index, i in enumerate(CNT):
			newatom = bgf.BgfAtom()
			newatom.x, newatom.y, newatom.z = i
			newatom.aNo = index
			newatom.aName = 'C' + str(index)
			newatom.rName = 'CNT'
			newatom.ffType = 'C_'
			newatom.chain = 'A'
			newatom.rNo = 1
			myBGF2.addAtom(newatom)

		for index, i in enumerate(WATER):
			newatom = bgf.BgfAtom()
			newatom.x, newatom.y, newatom.z = i
			newatom.aNo = index + n_CNT
			newatom.aName = 'O'
			newatom.rName = 'WAT'
			newatom.ffType = 'OW'
			newatom.chain = 'O'
			newatom.rNo = 100
			myBGF2.addAtom(newatom)

		myBGF2.REMARK.append('TIMESTEP ' + str(timestep))
		myBGF2.PERIOD = "111"
		myBGF2.AXES = "ZYX"
		myBGF2.SGNAME = "P 1                  1    1"
		myBGF2.CELLS = [-1, 1, -1, 1, -1, 1]
		myBGF2.CRYSTX = [boxsize[0], boxsize[1], boxsize[2], 90.0, 90.0, 90.0]
		myBGF2.saveBGF(temp_dir + bgf_file.split(".bgf")[0] + "." + str(timestep) + ".bgf")


		### move CM of CNT to origin
		CNT = CNT - boxsize / 2.0
		WATER = WATER - boxsize / 2.0


		### how the CNT is lying across the periodic box?
		for atom in CNT:
			min_x_CNT, min_y_CNT, min_z_CNT = CNT.min(axis=0)
			max_x_CNT, max_y_CNT, max_z_CNT = CNT.max(axis=0)
		
		margin = 5.0

		### copy water molecules
		max_x_axis = int(round((max_x_CNT+margin)/boxsize[0]))
		min_x_axis = int(round((min_x_CNT-margin)/boxsize[0]))
		max_y_axis = int(round((max_y_CNT+margin)/boxsize[1]))
		min_y_axis = int(round((min_y_CNT-margin)/boxsize[1]))
		max_z_axis = int(round((max_z_CNT+margin)/boxsize[2]))
		min_z_axis = int(round((min_z_CNT-margin)/boxsize[2]))

		replNum = (min_x_axis, max_x_axis, min_y_axis, max_y_axis, min_z_axis, max_z_axis)
		copyNum = 0

		remark = "";
		if watercopy:
			for x in range(min_x_axis, max_x_axis + 1):
				remark += "x: "
				if x != 0: 
					WATER2 = np.copy(WATERONLY)
					dx = np.array([[boxsize[0]*x, 0, 0]])
					WATER2 = WATER2 + dx
					WATER = np.concatenate((WATER, WATER2))
					copyNum += 1
					remark += str(x) + " "

			for y in range(min_y_axis, max_y_axis + 1):
				remark += "y: "
				if y != 0: 
					WATER2 = np.copy(WATERONLY)
					dy = np.array([[0, boxsize[1]*y, 0]])
					WATER2 = WATER2 + dy
					WATER = np.concatenate((WATER, WATER2))
					copyNum += 1
					remark += str(y) + " "

			for z in range(min_z_axis, max_z_axis + 1):
				remark += "z: "
				if z != 0: 
					WATER2 = np.copy(WATERONLY)
					dz = np.array([[0, 0, boxsize[2]*z]])
					WATER2 = WATER2 + dz
					WATER = np.concatenate((WATER, WATER2))
					copyNum += 1
					remark += str(z) + " "


		'''
		### WATER distribution check
		x_axis = np.linspace(WATER[:,0].min(), WATER[:,0].max(), 10)
		y_axis = np.linspace(WATER[:,1].min(), WATER[:,1].max(), 10)
		z_axis = np.linspace(WATER[:,2].min(), WATER[:,2].max(), 10)
		x_hist, _ = np.histogram(WATER[:,0], x_axis, normed=True)
		y_hist, _ = np.histogram(WATER[:,1], y_axis, normed=True)
		z_hist, _ = np.histogram(WATER[:,2], z_axis, normed=True)
		remark += " xdist: "
		for i in x_hist:
			remark += "{0:6.3f}".format(i)
		remark += " ydist: "
		for i in y_hist:
			remark += "{0:6.3f}".format(i)
		remark += " zdist: "
		for i in z_hist:
			remark += "{0:6.3f}".format(i)
		'''


		### MI of CNT calculation
		for atom in CNT:
			Ixx += (atom[1]**2 + atom[2]**2) / N_CNT
			Iyy += (atom[0]**2 + atom[2]**2) / N_CNT
			Izz += (atom[0]**2 + atom[1]**2) / N_CNT
			Ixy -= (atom[0] * atom[1]) / N_CNT
			Ixz -= (atom[0] * atom[2]) / N_CNT
			Iyz -= (atom[1] * atom[2]) / N_CNT
			
		I = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])	# the moment of inertia tensor
		eigval, eigvec = np.linalg.eig(I)	# eigval[0] is the minimum among the values.
		U = np.matrix(eigvec)
		Ut = U.T


		### box rotation
		for atom in CNT:
			v = np.matrix([atom[0], atom[1], atom[2]]).T
			Uv = Ut * v
			atom[0] = float(Uv[2]); atom[1] = float(Uv[1]); atom[2] = float(Uv[0])	# CNT rotation

		for atom in WATER:
			v = np.matrix([atom[0], atom[1], atom[2]]).T
			Uv = Ut * v
			atom[0] = float(Uv[2]); atom[1] = float(Uv[1]); atom[2] = float(Uv[0])	# water rotation


		### save coordinates to BGF after rotation
		myBGF2 = bgf.BgfFile()
		for index, i in enumerate(CNT):
			newatom = bgf.BgfAtom()
			newatom.x, newatom.y, newatom.z = i
			newatom.aNo = index
			newatom.aName = 'C' + str(index)
			newatom.rName = 'CNT'
			newatom.ffType = 'C_'
			newatom.chain = 'A'
			newatom.rNo = 1
			myBGF2.addAtom(newatom)

		for index, i in enumerate(WATER):
			newatom = bgf.BgfAtom()
			newatom.x, newatom.y, newatom.z = i
			newatom.aNo = index + n_CNT
			newatom.aName = 'O'
			newatom.rName = 'WAT'
			newatom.ffType = 'OW'
			newatom.chain = 'O'
			newatom.rNo = 100
			myBGF2.addAtom(newatom)

		myBGF2.REMARK.append('TIMESTEP ' + str(timestep))
		myBGF2.PERIOD = "111"
		myBGF2.AXES = "ZYX"
		myBGF2.SGNAME = "P 1                  1    1"
		myBGF2.CELLS = [-1, 1, -1, 1, -1, 1]
		myBGF2.CRYSTX = [boxsize[0], boxsize[1], boxsize[2], 90.0, 90.0, 90.0]
		myBGF2.saveBGF(temp_dir + bgf_file.split(".bgf")[0] + "." + str(timestep) + ".rot.bgf")


		### CNT height
		_, _, min_z_CNT = CNT.min(axis=0)
		_, _, max_z_CNT = CNT.max(axis=0)
		height_CNT = max_z_CNT - min_z_CNT


		### CNT radius
		x_CNT, y_CNT, z_CNT = np.mean(CNT, axis=0)
		x_std_CNT, y_std_CNT, z_std_CNT = np.std(CNT, axis=0)
		if x_std_CNT > 1.0 or y_std_CNT > 1.0:
			remark += ""


		### radius of CNT
		l_r_CNT = []; 
		for atom in CNT:
			l_r_CNT.append( math.sqrt( (atom[0] - x_CNT)**2 + (atom[1] - y_CNT)**2 ))

		r_CNT = np.mean(l_r_CNT)
		std_r_CNT = np.std(l_r_CNT)


		### get water molecules in CNT
		# inside the CNT := min_z_CNT <= z <= max_z_CNT and (x - (x_diff/2))**2 + (y - (y_diff/2))**2 < r_CNT**2
		# aNo_WAT_O_atoms: molecules which O atom is within CNT
		# we don't need to calculate H atoms. Let's consider only O atoms
		#####
		margin = 0.0;	# water molecules far from the margin will be only considered
		WAT_in_CNT = []
		for atom in WATER:
			dist_sq = (atom[0] - x_CNT)**2 + (atom[1] - y_CNT)**2
			if min_z_CNT + margin <= atom[2] and atom[2] <= max_z_CNT - margin and dist_sq < r_CNT**2:
				WAT_in_CNT.append(atom);

		n_WAT_in_CNT = len(WAT_in_CNT)


		'''
		### WATER bad contact check: copy-failure-proof: NEED CORRECTION
		WAT1 = []; WAT2 = [];
		for index1, i in enumerate(WAT_in_CNT):
			for j in WATER[index1:]:
				WAT1.append(i); WAT2.append(j);
		WAT1 = np.array(WAT1); WAT2 = np.array(WAT2)

		min_dists = np.min(np.dstack(((WAT1 - WAT2) % boxsize, (WAT2 - WAT1) % boxsize)), axis = 2)
		dists = np.sqrt(np.sum(min_dists ** 2, axis = 1))
		for d in dists:
			if d > 0 and d < 1.0:
				remark += "Bad contacts" + str(d)
				continue;
		'''



		d = "{0:8.3f}"; e = "{0:6.1f}"; output = "{0:<10}".format(timestep) + str(n_WAT_in_CNT) + ' | ' + d.format(r_CNT) + d.format(std_r_CNT) + ' | ' + e.format(boxsize[0]) + e.format(boxsize[1]) + e.format(boxsize[2]) + ' | ' + e.format(min_x_CNT) + e.format(max_x_CNT) + e.format(min_y_CNT) + e.format(max_y_CNT) + e.format(min_z_CNT) + e.format(max_z_CNT) + ' | ' + e.format(height_CNT) + ' | ' + str(replNum) + '\t' + str(copyNum) + '\t' + str(len(WATER)) + '\t' + str(remark) + '\n'
		ftemp.write(output)
		sys.stdout.flush()

		t2 = time.time()	# time mark
		elapsed_time = t2 - t1;
		processed_step += 1;


	print('')
	ftemp.close()
	print("Numbers of water molecules are written in " + out_file + " ..Done.")


	return 1

	### end of function


if __name__ == "__main__":
	bgf_file = ""; trj_file = ""; ff_file = ""; n_step = 0; watercopy = False;

	options, args = getopt.getopt(sys.argv[1:], 'hb:t:n:Co:f:', ['help', 'bgf=', 'trj=', 'step=', 'watercopy=', 'out=', 'ff='])

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
		elif option in ('-n', '--step'):
			n_step = int(value)
		elif option in ('-C', '--watercopy'):
			watercopy = True
		elif option in ('-o', '--out'):
			out_file = str(value)
		elif option in ('-f', '--ff'):
			ff_file = str(value)
	        elif option == NULL:
			print(usage)
			sys.exit(0)
	
	# more options: resname CNT, resname solvent, save BGF

	# main call
	countWater(bgf_file, trj_file, n_step, watercopy, ff_file)
