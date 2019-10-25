#!/home/noische/program/python27/bin/python
"""
crosslink.py
Original: Mar 14 2011 In Kim

Modifications:
- 111202
	adjacent crosslinking is prohibited. (adjacency > 8)
	naming atoms in cross-linker is changed.
	an H atom connected to head_C of the cross-linker will be randomly chosen and detached for connection.
	cross-linker will be aligned toward tail_N.
	the H atom in tail_xlnk which is the nearest from tail_N will be selected for connection.
	treating charges are changed.
- 110902
	distance criteria is adjusted to <> 5.0 6.0
- 110822
	removing bad contacts are disabled.
- 110810
	now calculates primary, secondary, and tertiary amines.

- 110630
	new filter added: adjacent nitrogens cannot be crosslinked. (criteria len(shortestpath) < 5)
- 110603
	createLammpsInput command is modified for "within-solvent crosslink" model
- 110602
	new patch method with patch_crosslink is applied with pair_coeff patches.
- 110518 
	Added the following to prevent hydrogen branch selection error
		if temp_list == []: continue;
- 110420
	MD procedure "min heat nvt" is changed to "min nvt".
"""

# Python Modules
import sys
import os
import string
import shutil
import random
import time
import copy
from types import *

# BGF Modules
sys.path.append("/home/noische/script")
sys.path.append("/home/noische/scripts")
from bgf import *
import bgftools
import bgfselect
import nutils as nu
from updateBGF import *
from centerBGF import *
from removeBadContacts import *
import patch_crosslink
import evaporatebgf

# Pizza.py Modules
sys.path.append("/home/noische/program/pizza-1Oct10")
from dump import *

# Cluster settings
from clusterSetting import *

# Globals
version = '111202'

def crosslink(bgf_file, ff_file, suffix, init_script_original_file, crosslink_script_original_file, t, ratio, out_file):

	# *** variables ***

	# used for calculating the number of 1", 2", and 3" amines
	l_init_amines = []; l_amines = [];

	# used for log file output
	output = ""

       	# initialization

	f_out_file = open(out_file, 'w')
	t1 = 0; t2 = 0;	# clock for measuring time for a cycle
	
	f_out_file.write("crosslink.py: version " + str(version) + "\n")
	f_out_file.write("" + "\n")
	f_out_file.write("Job started at " + time.asctime(time.gmtime()) + " on " + os.environ["HOSTNAME"] + " by " + os.environ["USER"] + "\n")
	f_out_file.write("Command executed at " + os.getcwd() + "\n")
	f_out_file.write("Requested options: " + str(sys.argv) + "\n")
	f_out_file.write("" + "\n")
	curr_dir = os.path.abspath(".")
	temp_dir = scratch_dir + "/" + suffix
	initial_bgf_file = bgf_file

	# LAMMPS executive command: defined in ~/scripts/clusterSettings.py
	lammps_parallel = mpi_command + " " + lammps_command

	n_crosslink = 0;
	f_out_file.flush()
	
	f_out_file.write("\n* Initialization" + "\n")

	# determine the initial numbers of primary, secondary, and tertiary amines
	l_init_amines = bgftools.getAmineGroupInfo(initial_bgf_file)
	f_out_file.write("\tNumber of primary, secondary, and tertiary amine groups: " + str(l_init_amines) + "\n")

	# scratch setting and file copy
	f_out_file.write("\tCreating scratch directory: " + temp_dir + "\n")
	if not os.path.isdir(temp_dir):
		os.makedirs(temp_dir)

	# managing hostfile: "COMMENTED"
	# Important Note: assume that "cat $PBS_NODEFILE > hostfile" is already executed in the .pbs file
	#hostfile = open(temp_dir + "/hostfile", 'w')
	#for i in range(0, 12): hostfile.write(os.environ["HOST"] + "\n")
	#hostfile.close()
	shutil.copy(initial_bgf_file, temp_dir)
	shutil.copy(ff_file, temp_dir)
	os.chdir(temp_dir)

	# LAMMPS input script files for crosslinking
	#init_script_original_file = "/home/noische/research/dendrimer/simulation/600/procedures/in.noinit"
	#init_script_original_file = "/home/noische/research/dendrimer/simulation/600/procedures/in.initialize"
	#crosslink_script_original_file = "/home/noische/research/dendrimer/simulation/600/procedures/in.xlnk_npt_1ps_k16_r28"
	#crosslink_script_original_file = "/home/noische/research/dendrimer/simulation/600/procedures/in.xlnk_npt_1ps_k8_r28"
	#crosslink_script_original_file = "/home/noische/research/dendrimer/simulation/600/procedures/in.xlnk_npt_1ps_k16_r35"

	shutil.copy(init_script_original_file, temp_dir)
	shutil.copy(crosslink_script_original_file, temp_dir)

	init_script_filename = init_script_original_file.split("/")[-1]	# initial script filename
	crosslink_script_filename = crosslink_script_original_file.split("/")[-1]	# crosslinking script filename

	createLammpsInput = "createLammpsInput.pl -b " + bgf_file + " -f " + ff_file + " -s " + suffix + " -o 'no shake' -t " + init_script_filename + " > /dev/null"
	os.system(createLammpsInput)
	initial_in_file = "in." + suffix
	
	# apply a patch for initialization
	patch_crosslink.pair_coeff_patch_wo_xlinker(initial_in_file, False)

	# first run
	t1 = time.time()
	runLammpsParallel = lammps_parallel + " -screen none -in " + initial_in_file + " -log " + os.path.join(temp_dir, suffix + ".log")
	f_out_file.write("\tRunning Lammps: see " + str(os.path.join(temp_dir, suffix + ".log")) + "\n")
	os.system(runLammpsParallel)
	t2 = time.time()
	f_out_file.write("\tElapsed time for the initialization: " + str(t2 - t1) + " sec\n")
	f_out_file.flush()
	
	# update coordinates
	initial_trj_file = suffix + "_init.lammpstrj"
	f_out_file.write("\tUpdating the BGF file " + initial_bgf_file + " with LAMMPS trajectory file to _trj_updated.bgf" + "\n")
	updatebgf(initial_bgf_file, initial_trj_file, "_trj_updated.bgf", True, silent=True);

	f_out_file.write("\tReassigning the origin of " + "_trj_updated.bgf" + " and saving " + "_centered.bgf" + "\n")
	#nu.shutup(); centerbgf("_trj_updated.bgf", "_centered.bgf", "", "com_center"); nu.say()
	nu.shutup(); os.system("centerBGF.pl -b _trj_updated.bgf -f " + ff_file + " -s _centered.bgf -c com_center"); nu.say()

	bgf_file = suffix + "1.bgf"
	shutil.copy("_centered.bgf", bgf_file)
	f_out_file.flush()
	
	# repeat
	if t == 0: t = 10;

	# crosslinking start!
	for i in range(1, t + 1):
		print("====== " + str(i) + " ======")
		t1 = time.time()
		N_prisecter = [0, 0, 0];	# N counter init.
		f_out_file.write("\n\n====== Cycle " + str(i) + " ======" + "\n")
	
		sname = suffix + str(i)
		bgf_file = sname + ".bgf"
		input_file = "in." + sname

		f_out_file.flush()

		# create input
		f_out_file.write("\n* Minimization, heating, and NVT dynamics for 1ps" + "\n")
		f_out_file.write("\tCreating Lammps input for " + bgf_file + "\n")
		################## IMPORTANT: min heat nvt -> min nvt for crosslinking ##################### 2011.4.20.

		# crosslinking with spring
		createLammpsInput = "createLammpsInput.pl -b " + bgf_file + " -f " + ff_file + " -s " + sname + " -o 'no shake' -t 'min " + crosslink_script_filename + "' > /dev/null"

		nu.shutup(); os.system(createLammpsInput); nu.say();

		##### apply a patch after the first step
		if n_crosslink != 0:
			patch_crosslink.pair_coeff_patch_w_xlinker("in." + sname, False)
		else:
			patch_crosslink.pair_coeff_patch_wo_xlinker("in." + sname, False)

		f_out_file.flush()

		# run lammps
		f_out_file.write("\tRunning LAMMPS simulation on " + sname + "\n")
		f_out_file.write("\tRunning Lammps: see " + str(os.path.join(temp_dir, sname + ".log")) + "\n")
		runLammpsParallel = lammps_parallel + " -screen none -in " + input_file + " -log " + sname + ".log"
		os.system(runLammpsParallel)
	
		# load trajectory
		f_out_file.write("\tLoading Trajectory." + "\n")
		initial_trj_file = sname + "_nvt.lammpstrj"
			
		# update bgf
		f_out_file.write("\tUpdating BGF file from trajectory." + "\n")
		#updatebgf(bgf_file, temp_trj_file, "_trj_updated.bgf", True, silent=True)
		updatebgf(bgf_file, initial_trj_file, "_trj_updated.bgf", True, silent=True)
	
		## align bgf periodic box	### added at 110719
		#bgftools.periodicMoleculeSort("_trj_updated.bgf", "_aligned.bgf", silent=True)

		# center bgf
		#nu.shutup(); centerbgf("_trj_updated.bgf", "_centered.bgf", "", "com_center"); nu.say()
		nu.shutup(); os.system("centerBGF.pl -b _trj_updated.bgf -f " + ff_file + " -s _centered.bgf -c com_center"); nu.say()
		shutil.copy("_centered.bgf", bgf_file)

		f_out_file.flush()

		# check the criteria from BGF file
		### 1. open the BGF file 
		f_out_file.write("\n* Crosslinking module execution" + "\n")

		clusterBGF = BgfFile(bgf_file)
		clusterBGF.renumber()
		f_out_file.write("\tChecking the connectivities in " + bgf_file + "\n")

		### 2. list up all N atoms
		N_index = [];		# REMARK: actually, this is aNo
		NN_dist_ok = [];
	
		for atom in clusterBGF.a:
			if string.strip(atom.ffType) == "N_3":
				N_index.append(atom.aNo)
				if "PRI" in atom.rName:
					N_prisecter[0] += 1
				elif "SEC" in atom.rName:
					N_prisecter[1] += 1
				elif "TER" in atom.rName:
					N_prisecter[2] += 1
		f_out_file.write("\tThe number of 1', 2', and 3' amine: " + str(N_prisecter) + "\n")
		N_number = len(N_index);

		criteria = 5.0;
	
		### 3. find N-N pairs which meet
		# 1st criteria: distance is between 5 and 6
		# 2nd criteria: only PRI or SEC amine
		# 3rd criteria: avoid repeating (removed)
		for x in range(N_number):
			for y in range(x + 1, N_number):
				if clusterBGF.distance(clusterBGF.a2i[N_index[x]], clusterBGF.a2i[N_index[y]]) < 6.0 \
					and clusterBGF.distance(clusterBGF.a2i[N_index[x]], clusterBGF.a2i[N_index[y]]) > 5.0 \
					and clusterBGF.getAtom(N_index[x]).rName != "TER" \
					and clusterBGF.getAtom(N_index[y]).rName != "TER": 
					NN_dist_ok.append([N_index[x], N_index[y]])

		### Additional Filtering (appended by 110628)
		##### This filter prevents the crosslinking of adjacent N-N pairs (self-pairing)

		"""
		#branchedoff_BGF = bgf.BgfFile()
		#branchedoff_BGF = copy.deepcopy(clusterBGF)
		#branchedoff_BGF = evaporatebgf.evaporate(branchedoff_BGF)
		#branchedoff_BGF = bgftools.getBackbone(branchedoff_BGF)
		#atompath = [];

		if NN_dist_ok != []:
			f_out_file.write("\tInitial crosslink candidate site: " + str(NN_dist_ok) + "\n")
			for pair in NN_dist_ok:
				f_out_file.write("\tChecking the path of " + str(pair))
				# if two candidate atom is in the same molecule, check the connectivity:
				#print(clusterBGF.getAtom(pair[0]))
				#print(clusterBGF.getAtom(pair[1]))
				pairaResName = int(clusterBGF.getAtom(pair[0]).rNo)
				pairbResName = int(clusterBGF.getAtom(pair[1]).rNo)

				if pairaResName == pairbResName:
					branchedoff_BGF = bgfselect.selectResNo(clusterBGF, pairaResName)
					atompath = [];
					try:
						atompath = bgftools.getShortestPath(branchedoff_BGF, pair[0], pair[1])
					except:
						atompath = [1, 2, 3, 4, 5, 6]	# it getShortestPath fails, then the two can be connected
						#nu.die("cannot calculate atom path from " + str(pair) + " at cycle " + str(i))
						#sys.exit(0)
	
					if not type(atompath) is NoneType:
						if len(atompath) < 5:
							f_out_file.write(".. An adjacent pair.. rejected!\n")
							NN_dist_ok.remove(pair)
						else:
							f_out_file.write(".. Same residue but far away.. accepted!\n")
				else:
					f_out_file.write(".. between separated molecule.. accepted!\n")
		"""
		f_out_file.write("\tCrosslinking candidate site: " + str(NN_dist_ok) + "\n" + "\tNow filtering out the candidate.\n")

		# now working on this 2011.11.30.inkim
		
		##### Filtering out the nearest adjacent and the second nearst adjacent cross-linking
		f_out_file.write("\tDeleting the nearest adjacent or second nearest adjacent from cross-link candidate site.\n")

		if NN_dist_ok != []:
			f_out_file.write("\tInitial crosslink candidate site: " + str(NN_dist_ok) + "\n")
			for pair in NN_dist_ok:
				branchedoff_BGF = bgf.BgfFile()	# initialize
				branchedoff_BGF = copy.deepcopy(clusterBGF)	# copy a cluster
				amine1 = branchedoff_BGF.getAtom(pair[0])
				amine2 = branchedoff_BGF.getAtom(pair[1])
				branchedoff_BGF = evaporatebgf.evaporate(branchedoff_BGF)	# delete water
				branchedoff_BGF = bgftools.getBackbone(branchedoff_BGF)		# delete hydrogens
				# REMARK: note that although branchedoff_BGF is edited, amine1 and amine2 points the atoms correctly.
				
				# for amine1 and amine2, count the third connected atom
				#          a1_left       a1_right
				# N --- C --- C <--- N ---> C --- C --- N
				#                  amine1
				amine1_conect = dict(); amine2_conect = dict()
				amine1_conect[amine1.aNo] = amine1.CONECT
				amine2_conect[amine2.aNo] = amine2.CONECT

				# 3-step connection for amine1
				for j in range(0, 2):
					d_temp = copy.deepcopy(amine1_conect)
					for ano in d_temp:
						for conect_ano in d_temp[ano]:
							tempatom = branchedoff_BGF.getAtom(conect_ano)
							amine1_conect[tempatom.aNo] = tempatom.CONECT
				
				# 3-step connection for amine2
				for j in range(0, 2):
					d_temp = copy.deepcopy(amine2_conect)
					for ano in d_temp:
						for conect_ano in d_temp[ano]:
							tempatom = branchedoff_BGF.getAtom(conect_ano)
							amine2_conect[tempatom.aNo] = tempatom.CONECT

				amine1_conect.update(amine2_conect)
				# now we get a small map for determining the connectivity
				amine1_to_amine2 = bgftools.findShortestPath(amine1_conect, amine1.aNo, amine2.aNo)
				f_out_file.write("\t\t* " + str(pair) + ": " + str(amine1_to_amine2) + "\n")
				if not type(amine1_to_amine2) is NoneType:
					if len(amine1_to_amine2) < 8:
						f_out_file.write("\t\tAdjacent pair.. will not be cross-linked.\n")
						NN_dist_ok.remove(pair)


		##### Filtering out the repeat
		# Tested with pair_candidate = [[1, 2], [1, 3], [1, 4], [3, 4], [3, 5], [5, 6], [7, 8], [9, 10], [9, 11]]
		if NN_dist_ok != []:
			result = [];
			pair_candidate = copy.deepcopy(NN_dist_ok)
			pair_candidate_flatten = cu.flatten(pair_candidate)
		
			for pair in pair_candidate:
				if pair_candidate_flatten.count(pair[0]) == 1 and pair_candidate_flatten.count(pair[1]) == 1:
					# unique pair in the candidate list
					result.append(pair)
					pair_candidate.remove(pair)	# unique pairs will be deleted from pair_candidate
		
			if pair_candidate != []:
				# grouping
				x = pair_candidate[0][0]; templist = []; grouped = [];
				for pair in pair_candidate:
					if pair[0] == x:
						templist.append(pair)
					else:
						grouped.append(templist)
						templist = []
						x = pair[0];
						templist.append(pair)
				grouped.append(templist)
			
				# shuffle pairs
				for group in grouped:
					random.shuffle(group)
				random.shuffle(grouped)
			
				# selecting a pair from group
				for group in grouped:
					pair = random.choice(group)
					result.append(pair)
			
				# removing repeats
				temp_repeated = []
				result_flatten = cu.flatten(result)
				for element in result_flatten:
					if result_flatten.count(element) != 1:	# if repeated..
						temp_repeated.append(element)
			
				for element in temp_repeated:
					templist = []
					for pair in result:
						if element in pair:
							templist.append(pair)
							result.remove(pair)
					if templist != []:
						selected = random.choice(templist)
						result.append(selected)
		
				NN_dist_ok = result


				##### end of filtering

		f_out_file.write("\t" + str(len(NN_dist_ok)) + " N-N pairs which meets the criteria ( " + str(criteria) + "A ) are found." + "\n")
		f_out_file.write("\t" + "The following N-N pairs will be crosslinked: " + str(NN_dist_ok) + "\n")
		f_out_file.flush()

		### 4. attach crosslinker
		# NN_dist_ok is the list that will be get crosslinked
		if len(NN_dist_ok) > 0:
			f_out_file.write("\n* Attaching the crosslinker to " + bgf_file + "\n")
			for pair in NN_dist_ok:
				f_out_file.write("\tLinking " + str(pair) + " Nitrogens: ")

				##### a) load a crosslinker
				xlinkerBGF_file = "/home/noische/research/dendrimer/xlinker/20111129/xlinker_new.bgf"
				#xlinkerBGF = BgfFile("/home/noische/research/dendrimer/structure/xlinker/xlnker.bgf")
				xlinkerBGF = BgfFile(xlinkerBGF_file)
				f_out_file.write("\tLoaded the crosslinker file " + str(xlinkerBGF_file) + "\n")
		
				head_xlnk = xlinkerBGF.a[0]	# head C of xlnk
				tail_xlnk = xlinkerBGF.a[2]	# tail C of xlnk
		
				##### b) select a hydrogen to be detached
				#    /- head_detach_H(M) /- head_xlnk_detach_H(H)      /- tail_detach_H(N)
				# head_N(O) --- head_xlnk(J) --- tail_xlnk(K) --- tail_N(P)
				#                                    \- tail_xlnk_detach_H(I)
				# Detach one H atom at random connected with J which is marked as X
				# Detach one H atom at random connected with K which is marked as Xmar
				head_N = clusterBGF.getAtom(pair[0])
				tail_N = clusterBGF.getAtom(pair[1])
				head_N.chain = "O"
				tail_N.chain = "P"
				f_out_file.write("(" + str(head_N.rName) + " and " + str(tail_N.rName) + " amines)\n")

				# change the crosslinker's rNo as same as head_N
				head_N_rNo = head_N.rNo
				for atom in xlinkerBGF.a:
					atom.rNo = head_N_rNo

				# choose which H atom will be detached from head_N
				temp_list = []
				for ano in head_N.CONECT:
					if clusterBGF.getAtom(ano).is_hydrogen():
						temp_list.append(ano)

				if temp_list == []: continue;

				head_H_detach_candidate_id = random.choice(temp_list)
				f_out_file.write("\tChoosing the hydrogen " + str(head_H_detach_candidate_id) + " from " + str(temp_list) + " attached to " + str(head_N.aNo) + "\n")

				# choose which H atom will be detached from tail_N
				temp_list = []
				for ano in tail_N.CONECT:
					if clusterBGF.getAtom(ano).is_hydrogen():
						temp_list.append(ano)

				if temp_list == []: continue;

				tail_H_detach_candidate_id = random.choice(temp_list)
				f_out_file.write("\tChoosing the hydrogen " + str(tail_H_detach_candidate_id) + " from " + str(temp_list) + " attached to " + str(tail_N.aNo) + "\n")

				head_detach_H = clusterBGF.getAtom(head_H_detach_candidate_id)	# H which will be detached
				tail_detach_H = clusterBGF.getAtom(tail_H_detach_candidate_id)	# also
				head_detach_H.chain = "M"
				tail_detach_H.chain = "N"
				head_detach_H_charge = head_detach_H.charge
				tail_detach_H_charge = tail_detach_H.charge
		
				# choose which hydrogen will be detached from J
				temp_list = [];
				for atom in xlinkerBGF.a:
					if atom.chain == "X":
						temp_list.append(atom.aNo)
				head_xlnk_detach_H_aNo = random.choice(temp_list);
				head_xlnk_detach_H = xlinkerBGF.getAtom(head_xlnk_detach_H_aNo)
				head_xlnk_detach_H.chain = "H"
				head_xlnk_detach_H_charge = head_xlnk_detach_H.charge
				f_out_file.write("\tH atom " + str(head_xlnk_detach_H_aNo) + " in the cross-linker will be deleted.\n")

				##### c) move xlinker near to the place: let N-H bond is 1.0 A: move C to 1.4 * N-H bond
				#delta = ( head_detach_H.x - head_xlnk.x, head_detach_H.y - head_xlnk.y, head_detach_H.z - head_xlnk.z )	# old method: moves head_xlnk to head_detach_H
				delta = ( 1.4 * head_detach_H.x - 0.4 * head_N.x - head_xlnk.x, 1.4 * head_detach_H.y - 0.4 * head_N.y - head_xlnk.y, 1.4 * head_detach_H.z - 0.4 * head_N.z - head_xlnk.z)	# new method
				moveBGF(xlinkerBGF, delta[0], delta[1], delta[2] )
		
				##### c2) rotate xlinker toward tail_N
				vec = (tail_N.x - head_N.x, tail_N.y - head_N.y, tail_N.z - head_N.z)
				xlinkerBGF = rotateBGF(xlinkerBGF, 1, 3, vec, 0)	# 0 means return as BgfFile

				##### c3) find which tail_detach_H is the nearst to tail_N
				temp_list = [];
				temp_dist = 100;
				for atom in xlinkerBGF.a:
					if atom.chain == "Y":
						temp_list.append(atom.aNo)
				tail_xlnk_detach_H_aNo = 0;

				for ano in temp_list:
					temp = clusterBGF.distance(clusterBGF.a2i[tail_N.aNo], xlinkerBGF.a2i[ano]);
					if temp < temp_dist:
						temp_dist = temp;
						tail_xlnk_detach_H_aNo = ano
				tail_xlnk_detach_H = xlinkerBGF.getAtom(tail_xlnk_detach_H_aNo)
				tail_xlnk_detach_H.chain = "I"
				tail_xlnk_detach_H_charge = tail_xlnk_detach_H.charge

				##### d) merge
				clusterBGF = clusterBGF.merge(xlinkerBGF, True)
				clusterBGF.saveBGF("_xlnk_temp.bgf")
				clusterBGF = BgfFile("_xlnk_temp.bgf")
		
				##### e) detach head_detach_H and connect head_N-head_xlnk and tail_N-tail_xlnk

				for atom in clusterBGF.a:
					if atom.chain == "O":
						head_N = atom
					if atom.rName == "XLK":
						if atom.chain == "J":
							head_xlnk = atom
					if atom.chain == "P":
						tail_N = atom
					if atom.rName == "XLK":
						if atom.chain == "K":
							tail_xlnk = atom
					if atom.chain == "H":
						head_xlnk_detach_H = atom
					if atom.chain == "I":
						tail_xlnk_detach_H = atom

				# disconnect Hydrogens from head_N and tail_N
				f_out_file.write("\tDisconnecting Hydrogen(aNo: " + str(head_detach_H.aNo) + ") from Nitrogen(aNo: " + str(head_N.aNo) + ") - Head" + "\n")
				f_out_file.write("\tDisconnecting Hydrogen(aNo: " + str(tail_detach_H.aNo) + ") from Nitrogen(aNo: " + str(tail_N.aNo) + ") - Tail" + "\n")
				clusterBGF.disconnect(clusterBGF.a2i[head_N.aNo], clusterBGF.a2i[head_detach_H.aNo])
				clusterBGF.disconnect(clusterBGF.a2i[tail_N.aNo], clusterBGF.a2i[tail_detach_H.aNo])

				# detach Hydrogens from head_N and tail_N
				if head_detach_H.aNo > tail_detach_H.aNo:
					f_out_file.write("\tDetaching Hydrogen atom: " + str(head_detach_H.aNo) + "\n")
					clusterBGF.delAtom(clusterBGF.a2i[head_detach_H.aNo])
					f_out_file.write("\tDetaching Hydrogen atom: " + str(tail_detach_H.aNo) + "\n")
					clusterBGF.delAtom(clusterBGF.a2i[tail_detach_H.aNo])
				else:
					f_out_file.write("\tDetaching Hydrogen atom: " + str(tail_detach_H.aNo) + "\n")
					clusterBGF.delAtom(clusterBGF.a2i[tail_detach_H.aNo])
					f_out_file.write("\tDetaching Hydrogen atom: " + str(head_detach_H.aNo) + "\n")
					clusterBGF.delAtom(clusterBGF.a2i[head_detach_H.aNo])

				# disconnect H-J and delete H 
				clusterBGF.disconnect(clusterBGF.a2i[head_xlnk.aNo], clusterBGF.a2i[head_xlnk_detach_H.aNo])
				clusterBGF.delAtom(clusterBGF.a2i[head_xlnk_detach_H.aNo])

				# disconnect K-I and delete I
				clusterBGF.disconnect(clusterBGF.a2i[tail_xlnk.aNo], clusterBGF.a2i[tail_xlnk_detach_H.aNo])
				clusterBGF.delAtom(clusterBGF.a2i[tail_xlnk_detach_H.aNo])

				# connect O-J
				clusterBGF.connect(clusterBGF.a2i[head_N.aNo], clusterBGF.a2i[head_xlnk.aNo])	# O-J connection

				# connect K-P
				clusterBGF.connect(clusterBGF.a2i[tail_N.aNo], clusterBGF.a2i[tail_xlnk.aNo])	# K-P connection

				# residue update
				f_out_file.write("\tUpdating residue information." + "\n")
				if head_N.rName == "PRI":
					head_N.rName = "SEC"
				elif head_N.rName == "SEC":
					head_N.rName = "TER"
				if tail_N.rName == "PRI":
					tail_N.rName = "SEC"
				elif tail_N.rName =="SEC":
					tail_N.rName = "TER"

				# chain update
				f_out_file.write("\tUpdating chain information." + "\n")
				head_N.chain = "V"
				tail_N.chain = "V"
				#head_xlnk.chain = "A"
				#tail_xlnk.chain = "A"
	
				# charge
				f_out_file.write("\tCalculating charges..\n")
				f_out_file.write("\tCharge of the amine group 1: " + str("{0:8.5f}".format(head_N.charge)) + "\n")
				f_out_file.write("\tCharge of the amine group 2: " + str("{0:8.5f}".format(tail_N.charge)) + "\n")
				f_out_file.write("\tCharge of the head C in the cross-link: " + str("{0:8.5f}".format(head_xlnk.charge)) + "\n")
				f_out_file.write("\tCharge of the tail C in the cross-link: " + str("{0:8.5f}".format(tail_xlnk.charge)) + "\n")
				f_out_file.write("\tCharge of the H detached from amine group 1: " + str("{0:8.5f}".format(head_detach_H_charge)) + "\n")
				f_out_file.write("\tCharge of the H detached from amine group 2: " + str("{0:8.5f}".format(tail_detach_H_charge)) + "\n")
				f_out_file.write("\tCharge of the H detached from head C in the cross-link: " + str("{0:8.5f}".format(head_xlnk_detach_H_charge)) + "\n")
				f_out_file.write("\tCharge of the H detached from tail C in the cross-link: " + str("{0:8.5f}".format(tail_xlnk_detach_H_charge)) + "\n")
				#head_N.charge += head_detach_H_charge
				#tail_N.charge += tail_detach_H_charge
				#tail_xlnk.charge += tail_xlnk_detach_H_charge
				#head_xlnk.charge += head_xlnk_detach_H_charge
				head_xlnk.charge += head_detach_H_charge
				tail_xlnk.charge += tail_detach_H_charge
				head_N.charge += head_xlnk_detach_H_charge
				tail_N.charge += tail_xlnk_detach_H_charge
				
				#leftoverCharge = clusterBGF.charge() / 2.0
				#head_N.charge -= leftoverCharge
				#tail_N.charge -= leftoverCharge
				leftoverCharge = clusterBGF.charge() 
				f_out_file.write("\tExcess Charge: " + str("{0:8.5f}".format(leftoverCharge)) + "\n")
				f_out_file.write("\tExcess Charge will be added to head N atom.\n")
				head_N.charge -= leftoverCharge
				
				f_out_file.write("\tModified charge of the amine group 1: " + str("{0:8.5f}".format(head_N.charge)) + "\n")
				f_out_file.write("\tModified charge of the amine group 2: " + str("{0:8.5f}".format(tail_N.charge)) + "\n")
				f_out_file.write("\tThe charge of the system is " + str("{0:8.5f}".format(clusterBGF.charge())) + "\n")
				f_out_file.write("" + "\n")

				n_crosslink += 1
				f_out_file.flush()
				##### The end of the crosslinking step

			##### Minimization
			clusterBGF.renumber()

			f_out_file.write("\tSaving the crosslinked file to " + "_linked.bgf" + "\n")
			clusterBGF.saveBGF("_linked.bgf")

			f_out_file.write("\n* Crosslinking step is finished. Entering minimizing step." + "\n")

			# remove bad contacts
			#f_out_file.write("\tRemoving bad contacts after crosslinking." + "\n")
			#nu.shutup(); removebadcontacts("_linked.bgf", "_removed.bgf", 2.0); nu.say();
		
			# create input for minimization
			#createLammpsInput = "createLammpsInput.pl -b _removed.bgf -f " + ff_file + " -s minimization -o 'no shake' -t 'min' > /dev/null"
			createLammpsInput = "createLammpsInput.pl -b _linked.bgf -f " + ff_file + " -s minimization -o 'no shake' -t 'min' > /dev/null"
			runLammpsParallel = lammps_parallel + " -screen none -in in.minimization -log minimization.log"
			f_out_file.write("\tMaking LAMMPS input script for minimization." + "\n")
			nu.shutup(); os.system(createLammpsInput); nu.say();

			# patch for minimization ## added at 110524
			if n_crosslink != 0:
				patch_crosslink.pair_coeff_patch_w_xlinker("in.minimization", False)
			else:
				patch_crosslink.pair_coeff_patch_wo_xlinker("in.minimization", False)
			f_out_file.write("\tApplied a patch on in.minimization file.\n")
		
			# run lammps
			f_out_file.write("\tMinimizing the structure." + "\n")
			f_out_file.write("\tRunning Lammps: see " + str(os.path.join(temp_dir, "minimization.log")) + "\n")
			os.system(runLammpsParallel)
		
			# load trajectory
			f_out_file.write("\tLoading Trajectory from the minimization result." + "\n")
			update_trj_file = "minimization_min.lammpstrj"
				
			# update bgf
			f_out_file.write("\tUpdating BGF file from the minimization result." + "\n")
			#updatebgf("_removed.bgf", "_temp.lammpstrj", "_trj_updated.bgf", silent=True)
			#updatebgf("_linked.bgf", "_temp.lammpstrj", "_trj_updated.bgf", True, silent=True)
			updatebgf("_linked.bgf", "minimization_min.lammpstrj", "_trj_updated.bgf", True, silent=True)
		
			# align bgf periodic box	### added at 110719
			bgftools.periodicMoleculeSort("_trj_updated.bgf", "_trj_updated.bgf", silent=True)

			# center bgf
			f_out_file.write("\tAligning BGF file." + "\n")
			#nu.shutup(); centerbgf("_trj_updated.bgf", "_centered.bgf", "", "com_center"); nu.say()
			nu.shutup(); os.system("centerBGF.pl -b _trj_updated.bgf -f " + ff_file + " -s _centered.bgf -c com_center"); nu.say()
			shutil.copy("_centered.bgf", suffix + str(i+1) + ".bgf")

			f_out_file.write("" + "\n")
			f_out_file.write("--- " + str(len(NN_dist_ok)) + " times of crosslinking occured at the cycle " + str(i) + "\n")
			f_out_file.write("--- End of the crosslinking process for cycle " + str(i) + "\n")
			f_out_file.write("--- See " + suffix + str(i+1) + ".bgf for output. \n" + "\n")
		else:
			f_out_file.write("--- Passing the crosslinking step in cycle " + str(i) + "\n" + "\n")
			shutil.copy(bgf_file, suffix + str(i+1) + ".bgf")
			f_out_file.flush()

		t2 = time.time()
		f_out_file.write("Elapsed time for the cycle: " + str(t2 - t1) + " sec\n")
		
		##### The end of the crosslinking cycle (1ps)

	f_out_file.write("\n\nJob completed at " + time.asctime(time.gmtime()) + "\n")
	f_out_file.write("Output structure on " + suffix + "_out.bgf" + "\n")
	f_out_file.write("The number of crosslinking through this job: " + str(n_crosslink) + "\n")
	shutil.copy(suffix + str(i+1) + ".bgf", suffix + "_out.bgf")
	os.system("cp -rp " + temp_dir + " " + curr_dir)

	##### End of crosslink.py


def rotateBGF(bgf_file, atom1, atom2, vec, out_file, silent=False):
	"""
rotateBGF: rotate v21 to the given vector (vec)
	"""
 
	# open
	if isinstance(bgf_file, bgf.BgfFile):
		myBGF = bgf_file
	else:
		if not silent: print("reading " + bgf_file + " ..")
		myBGF = bgf.BgfFile(bgf_file)

	a1 = myBGF.a[myBGF.a2i[atom1]];
	a2 = myBGF.a[myBGF.a2i[atom2]];

	v1 = (a1.x, a1.y, a1.z)

	## move atom1 to origin
	#for atom in myBGF.a:
	#	atom.x -= v1[0]
	#	atom.y -= v1[1]
	#	atom.z -= v1[2]

	v21 = [ (a2.x - a1.x), (a2.y - a1.y), (a2.z - a1.z) ]

	# rotate v21 to x axis
	u1 = v21 / np.linalg.norm(v21)	##

	v2 = [u1[1], -u1[0], 0]
	u2 = v2 / np.linalg.norm(v2)	##

	v3 = np.cross(u1, u2)
	u3 = v3 / np.linalg.norm(v3)	##

	U = np.array([[u1[0], u1[1], u1[2]], [u2[0], u2[1], u2[2]], [u3[0], u3[1], u3[2]]])

	# V = vec -> (1, 0, 0)
	# invV = inverse(V)
	V1 = vec / np.linalg.norm(vec)
	vec2 = [V1[1], -V1[0], 0]
	V2 = vec2 / np.linalg.norm(vec2)
	vec3 = np.cross(V1, V2)
	V3 = vec3 / np.linalg.norm(vec3)
	V = np.array([[V1[0], V1[1], V1[2]], [V2[0], V2[1], V2[2]], [V3[0], V3[1], V3[2]]])
	invV = np.linalg.inv(V)

	# rotate all atoms
	for atom in myBGF.a:
		a = np.matrix([atom.x, atom.y, atom.z]).T
		b = U*a
		c = invV*b
		atom.x = float(c[0])
		atom.y = float(c[1])
		atom.z = float(c[2])

	## save
	if isinstance(out_file, str):
		if not silent: print("saving information to " + out_file + " ..")
		myBGF.saveBGF(out_file)
		return 1;
	else:
		return myBGF;

	### end of function


if __name__ == '__main__':

	option = ""; args = ""; bgf_file = ""; ff_file = ""; out_file = "";
	suffix = ""; nodes = ""; t = 0; ratio = 0.0;
	init_script_original_file = ""; crosslink_script_original_file = "";
	usage = """
Usage: crosslink.py -b bgfFile -f forcefield -s suffix -c xlnk_dist_criteria -t n_cycles -o log_file -r ratio
	"""

	if len(sys.argv) < 2:
		print(usage);
		sys.exit(0)

	options, args = getopt.getopt(sys.argv[1:], 'hb:f:s:i:c:t:o:r:', ['help','bgf=','forcefield=','suffix=','input=','crosslink=','time=','output=','ratio='])
	for option, value in options:
		if option in ('-h', '--help'):
			print usage; sys.exit(0)
		elif option in ('-b', '--bgf'):
			bgf_file = value
		elif option in ('-f', '--forcefield'):
			ff_file = value
		elif option in ('-s', '--suffix'):
			suffix = value
		elif option in ('-i', '--input'):
			init_script_original_file = value
		elif option in ('-c', '--control'):
			crosslink_script_original_file = value
		elif option in ('-t', '--time'):
			t = int(value)
		elif option in ('-o', '--option'):
			out_file = value
		elif option in ('-r', '--ratio'):
			ratio = float(value)
		elif option in (''):
			print usage; sys.exit(0)

		if suffix == "": suffix = "crosslink"
		if out_file == "": out_file = suffix + "_log.out"
		if ff_file == "": ff_file = "/home/noische/ff/dreiding-den.par"
		#if ratio == 0.0: ratio == 0.5

	crosslink(bgf_file, ff_file, suffix, init_script_original_file, crosslink_script_original_file, t, ratio, out_file)
