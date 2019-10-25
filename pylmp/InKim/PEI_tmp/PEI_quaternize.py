#!/home/noische/program/python27/bin/python
"""
quaternize.py
Original: Apr 19 2011 In Kim
"""

# Python Modules
import sys
import os
import string
import random
import getopt
import time

# Custom Modules
import bgf
import nutils as nu
import bgftools

# Globals
version = '120728'

def protonate(bgf_file, ff_file, probability, out_file, log_file, silent=False):
	"""
protonate(bgf_file, ff_file, suffix, nodes, criteria, t, out_file):
	Protonate (performs alkylation of primary amines) the given bgf file.

	To-do:
		clean temp file management
		display appropriate messages
		server environments
	"""

	# initialization

	# open bgffile
	myBGF = bgf.BgfFile(bgf_file)

	# find all tertiary amines
	N_index = [];
	last_hetatm_aNo = 0;
	for atom in myBGF.a:
		if string.strip(atom.ffType) == "N_3":
			if "PRI" in atom.rName:
				N_index.append(atom.aNo)

	N_number = len(N_index);	# the number of all tertiary amines
	if not silent: print(str(N_number) + " primary amines are found.")

	# find the last HETATM (insertion position)
	for atom in myBGF.a:
		if atom.aTag != 1:
			last_hetatm_aNo = atom.aNo
			break;
	print(last_hetatm_aNo)

	# choose the primary amines randomly according to the probability
	if not probability == 1.0:
		N_choose = int(N_number * probability)
		N_index = random.sample(N_index, N_choose)

	if not silent: print(str(len(N_index)) + " primary amines will be protonated.")

	# add an hydrogen
	if len(N_index) > 0:
		for index, aNo in enumerate(N_index):

			if not silent: print("*** Primary amine " + str(aNo) + " ***")
			atom_Npri = myBGF.getAtom(aNo)	# N from the N_index list

			# new proton on primary amine
			proton = bgf.BgfAtom();
			proton.x = atom_Npri.x + 1;
			proton.y = atom_Npri.y + 1;
			proton.z = atom_Npri.z + 1;
			proton.aName = "H+_" + str(index);
			proton.rName = "PRI";
			proton.chain = "A";
			proton.aTag = 1;
			proton.rNo = atom_Npri.rNo;
			proton.ffType = "H_";
			proton.aNo = last_hetatm_aNo

			myBGF.addAtom(proton, last_hetatm_aNo-1)
			proton = myBGF.getAtom(last_hetatm_aNo)		# since (proton) != (proton add in myBGF)
			#proton.aName = "H++"
			myBGF.connect(myBGF.a2i[atom_Npri.aNo], myBGF.a2i[proton.aNo])		# connect

			# charge change: NH3 - CH2 - CH2 - .....
			# 1. NH3 part
			N_CONECT_aNo = atom_Npri.CONECT;
			N_CONECT_C_aNo = [];
			N_CONECT_H_aNo = [];
			for ano in N_CONECT_aNo:
				atom = myBGF.getAtom(ano)
				if "C" in atom.ffType:
					N_CONECT_C_aNo.append(atom.aNo)
				elif "H" in atom.ffType:
					N_CONECT_H_aNo.append(atom.aNo)
				else:
					nu.die(str(atom))

			if len(N_CONECT_C_aNo) != 1:
				nu.die("Number of C atoms connected with primary amine " + str(atom_Npri.aNo) + " mismatch!")
			if len(N_CONECT_H_aNo) != 3:
				nu.die("Number of H atoms connected with primary amine " + str(atom_Npri.aNo) + " mismatch!")

			atom_C2 = myBGF.getAtom(N_CONECT_C_aNo[0])
			atom_HN1 = myBGF.getAtom(N_CONECT_H_aNo[0])
			atom_HN2 = myBGF.getAtom(N_CONECT_H_aNo[1])
			atom_HN3 = myBGF.getAtom(N_CONECT_H_aNo[2])

			# update charges
			old_NH3_charge = atom_Npri.charge + atom_HN1.charge + atom_HN2.charge + atom_HN3.charge

			atom_Npri.charge = -0.72540
			atom_HN1.charge = 0.43330;
			atom_HN2.charge = 0.43330;
			atom_HN3.charge = 0.43330;
			new_NH3_charge = atom_Npri.charge + atom_HN1.charge + atom_HN2.charge + atom_HN3.charge
			NH3_charge_change = new_NH3_charge - old_NH3_charge;
			if not silent: print("NH3 site charge change: " + str(atom_Npri.aNo) + " : " + str(NH3_charge_change))


			# 2. first CH2 part
			C2_CONECT_aNo = atom_C2.CONECT;
			C2_CONECT_C_aNo = [];
			C2_CONECT_H_aNo = [];
			for ano in C2_CONECT_aNo:
				atom = myBGF.getAtom(ano)
				if "H" in atom.ffType:
					C2_CONECT_H_aNo.append(atom.aNo)
				elif "C" in atom.ffType and ano != atom_C2.aNo:
					C2_CONECT_C_aNo.append(atom.aNo)

			if len(C2_CONECT_H_aNo) != 2:
				nu.die("Number of H atoms connected with primary amine " + str(atom_C2.aNo) + " mismatch!")
			if len(C2_CONECT_C_aNo) != 1:
				nu.die("Number of C atoms connected with primary amine " + str(atom_C2.aNo) + " mismatch!")
				
			atom_C3 = myBGF.getAtom(C2_CONECT_C_aNo[0])
			atom_HC21 = myBGF.getAtom(C2_CONECT_H_aNo[0])
			atom_HC22 = myBGF.getAtom(C2_CONECT_H_aNo[1])

			old_CH2_charge = atom_C2.charge + atom_HC21.charge + atom_HC22.charge
			#if not silent: print((atom_C2.charge, atom_HC21.charge, atom_HC22.charge))
			#if not silent: print("old CH2 site charge change: " + str(atom_Npri.aNo) + " : " + str(old_CH2_charge))

			atom_C2.charge = -0.20790
			atom_HC21.charge = 0.24090
			atom_HC22.charge = 0.24090
			new_CH2_charge = atom_C2.charge + atom_HC21.charge + atom_HC22.charge
			#if not silent: print("new CH2 site charge change: " + str(atom_Npri.aNo) + " : " + str(new_CH2_charge))
			CH2_charge_change = new_CH2_charge - old_CH2_charge
			if not silent: print("CH2 site charge change: " + str(atom_Npri.aNo) + " : " + str(CH2_charge_change))

			delta_charge = NH3_charge_change + CH2_charge_change
			if not silent: print("Primary amine site charge change: " + str(atom_Npri.aNo) + " : " + str(delta_charge))

			# 3. second CH2 part
			old_C3_charge = atom_C3.charge
			if not silent: print("Old C3 charge: " + str(atom_C3.charge))
			atom_C3.charge = atom_C3.charge + (1 - delta_charge)
			if not silent: print("New C3 charge: " + str(atom_C3.charge))
			if not silent: print("Total charge: " + " : " + str(myBGF.charge()))


		# renumber
		myBGF.renumber()

		# save
		myBGF.REMARK.insert(0, "Protonated by " + os.path.basename(sys.argv[0]) + " by " + os.environ["USER"] + " on " + time.asctime(time.gmtime()))
		myBGF.saveBGF(bgf_file[:-4] + "_protonated.bgf")
		if not silent: print("saving file " + bgf_file[:-4] + "_protonated.bgf")

	else:
		if not silent: print("No amine groups to be protonated"); sys.exit(0);

	if not silent: print("Done.")
	# minimize the structure
	#createLammpsInput = "~tpascal/scripts/createLammpsInput.pl -b _prot.bgf -f " + ff_file + " -s minimization -o 'no shake' -t 'min' > /dev/null"
	#lammps_parallel = "/home/caltech/openmpi-1.3.3/bin/mpirun -np 1 /home/noische/program/lammps/lammps-14Jan11/src/lmp_qch"
	#lammps_parallel = "/opt/mpi/intel/mpich2-1.4/bin/mpirun -np 16 -f hostfile /opt/applic/lammps/bin/lmp_kdft" 
	#runLammpsParallel = lammps_parallel + " -screen none -in in.minimization -log minimization.log"
	
	#nu.shutup(); 
	#os.system(createLammpsInput); 
	#os.system(runLammpsParallel); 

	# load the minimized trajectory
	#update_trj_file = "minimization_min.lammpstrj"
	#temp_trj_file = "_temp.lammpstrj"
	#dumpinfo = dump(update_trj_file)
	#l_timestep = dumpinfo.time()
	#last_timestep = l_timestep[-1]
	#dumpinfo.tselect.one(last_timestep)
	#dumpinfo.unwrap()
	#dumpinfo.write("_temp.lammpstrj")

	# update bgf with the trajectory
	#updatebgf("_prot.bgf", "_temp.lammpstrj", out_file);
	#nu.say();
	
	##### end of the quaternization #####


if __name__ == '__main__':

	option = ""; args = ""; bgf_file = ""; ff_file = ""; probability = 0; out_file = ""; log_file = "";
	usage = """
Usage: quaternize.py -b bgfFile -f forcefield -p probability -o outfile 
	"""

	if len(sys.argv) < 2:
		print(usage)
		sys.exit(0)

	options, args = getopt.getopt(sys.argv[1:], 'hb:f:p:o:', ['help','bgf=','forcefield=','probability=','output='])
	for option, value in options:
		if option in ('-h', '--help'):
			print usage; sys.exit(0)
		elif option in ('-b', '--bgf'):
			bgf_file = value
		elif option in ('-f', '--forcefield'):
			ff_file = value
		elif option in ('-p', '--probability'):
			probability = float(value)
		elif option in ('-o', '--option'):
			out_file = value
		elif option in ('-l', '--log'):
			log_file = value
		elif option in (''):
			print usage; sys.exit(0)

		if probability == 0: probability = 1
		if out_file == "": out_file = bgf_file[:-4] + "_protonated.bgf"

	protonate(bgf_file, ff_file, probability, out_file, log_file, False)
