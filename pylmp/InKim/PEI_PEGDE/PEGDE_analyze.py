#!/home/noische/python

import sys
import os
import string
import copy

import bgf
import bgftools
import nutils as nu
import dreiding
import numpy
from itertools import dropwhile

usage = """
PEGDE_analyze.py bgffile qcoutfile
"""

#getcontext().prec = 7	# for charge summation in bgf in decimal package; deprecated.


def rindex(lst, item):
    def index_ne(x):
        return lst[x] != item
    try:
        return dropwhile(index_ne, reversed(xrange(len(lst)))).next()
    except StopIteration:
        raise ValueError, "rindex(lst, item): item not in list"


def read_data_from_qcout(qcout):

	# read
	f = open(qcout)
	fsp = f.readlines()

	# parse method and jobtype;;;; later

	# find the last shot
	cycles = [i for i, l in enumerate(fsp) if 'Optimization Cycle' in l]

	if len(cycles) == 0:
		nu.die("Please check the Q-Chem output file: is this an opt job output?")

	fsp = fsp[cycles[-1]:]	# cut the last cycle of the output

	# atom types and coordinates
	x = [i for i, l in enumerate(fsp) if 'Coordinates' in l][0]
	y = [i for i, l in enumerate(fsp) if 'Number of degrees of freedom' in l][0]
	
	coord = fsp[x+2:y]	# slice coordinates
	for i, l in enumerate(coord):
		coord[i] = coord[i].split()	# trim

	# charge
	x = [i for i, l in enumerate(fsp) if 'chelpg' in l.lower()]

	useMulliken = False;
	if len(x) == 0:
		print("** No ChElPG. Using Mulliken for charges.")
		useMulliken = True;
	else:
		print("** Using ESP charges from ChElPG calculation.")

	if useMulliken:
		x = [i for i, l in enumerate(fsp) if 'Net Atomic Charges' in l][0]
		y = [i for i, l in enumerate(fsp) if 'Sum of atomic charges' in l][0]
	else:
		x = [i for i, l in enumerate(fsp) if 'ChElPG Net Atomic Charges' in l][0]
		y = [i for i, l in enumerate(fsp) if 'Sum of atomic charges' in l][1]	# if there is ChElPG, then y returns two values

	charge = fsp[x+4:y-1]
	for i, l in enumerate(charge):
		charge[i] = charge[i].split()

	if len(charge) != len(coord):
		nu.warn("Error on reading qcout file. Please check the number of atoms.")
	
	return coord, charge; 


def calculate_mi_rotate(myBGF):
	"""
	require myBGF as bgf object and returns rotated bgf object
	get from "/qcfs/noische/scripts/CNT_alignCNT.py"
	"""
	# initialize for the moment of inertia and the center of mass calculation
	U = 0; Ut = 0; Uv = 0;
	Ixx = 0.0; Ixy = 0.0; Ixz = 0.0;
	Iyx = 0.0; Iyy = 0.0; Iyz = 0.0;
	Izx = 0.0; Izy = 0.0; Izz = 0.0;
	Mx = 0.0; My = 0.0; Mz = 0.0; m = 0.0;

	# load mass from DREIDING2
	FF = dreiding.loadFF("")
	atominfo = dreiding.loadAtomTypes(FF)
	
	# calculate the moment of inertia (MI) and the center of mass (COM) of CNT from "myBGF"
	for atom in myBGF.a:
		fftype_key = string.strip(atom.ffType)
		mi = float(atominfo[fftype_key]['MASS'])
		Mx += atom.x * mi
		My += atom.y * mi
		Mz += atom.z * mi
		Ixx += (atom.y**2 + atom.z**2) * mi
		Iyy += (atom.x**2 + atom.z**2) * mi
		Izz += (atom.x**2 + atom.y**2) * mi
		Ixy -= (atom.x * atom.y) * mi
		Ixz -= (atom.x * atom.z) * mi
		Iyz -= (atom.y * atom.z) * mi
		m += mi 

	Mx /= m; My /= m; Mz /= m;
	I = numpy.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])	# the moment of inertia tensor

	# transpose for "all atoms in BGF": move COM as origin
	for atom in myBGF.a:
		atom.x -= Mx
		atom.y -= My
		atom.z -= Mz

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

	return myBGF;


def main():
	# remark: play with /qcfs/noische/research/PEI/201410/EP-PEG-EP/EP-PEG4-EP/result/EP-PEG4-EP_opt_b3lyp_6-311pgss.out

	if len(sys.argv) < 2:
		nu.die("Please specify relevant files. Usage: " + usage)
	bgffile = sys.argv[1]
	qcoutfile = sys.argv[2]

	# get coord and charge
	print("** fetching coordinate and charge from the qcout file.")
	coord, charge = read_data_from_qcout(qcoutfile)

	# update coord and charge if name and aNo matches
	# assume that bgf is well converted from something....;;
	print("** applying coordinate and charge from the qcout file.")
	bgf1 = bgf.BgfFile(bgffile)
	for i, atom in enumerate(bgf1.a):
		if atom.aNo == int(coord[i][0]) and atom.aNo == int(charge[i][0]):
			atom.x = float(coord[i][2])
			atom.y = float(coord[i][3])
			atom.z = float(coord[i][4])
			atom.charge = float(charge[i][2])


	### calculate MI and align molecule
	bgf1 = calculate_mi_rotate(bgf1)


	### find epoxide group
	aNo_O = []; aNo_epoxide = []; aNo_O_epoxide = []; aNo_identified = []; aNo_ethylenes = []
	for atom in bgf1.a:
		if "O" in atom.aName:
			aNo_O.append(atom.aNo)

	# check C-O-C arm connectivity
	for i, ano in enumerate(aNo_O):
		atomO = bgf1.getAtom(aNo_O[i])
		atomC1 = bgf1.getAtom(atomO.CONECT[0])
		atomC2 = bgf1.getAtom(atomO.CONECT[1])
		if bgf.is_bonded(atomC1, atomC2):
			aNo_epoxide.append([ano, atomO.CONECT[0], atomO.CONECT[1]])	# epoxide O and two C atoms

	if len(aNo_epoxide) != 2:
		nu.warn("More than two epoxide groups are found.")

	#!print(aNo_epoxide)

	# mark for head: the first element will be the head
	aNo_head_O = aNo_epoxide[0][0]
	aNo_start = 0;
	aNo_end = 0;

	# find the nearst O from the head
	for ano in aNo_O:
		d = len(bgftools.getShortestPath(bgf1, aNo_head_O, ano))
		if d > 1 and d < 5:
			aNo_epoxide[0].append(ano)
			aNo_start = ano	# starting point for ethylene group

	# the C between head O and second O
	for atom in bgf1.a:
		d = len(bgftools.getShortestPath(bgf1, aNo_head_O, atom.aNo))
		if d == 3 and "C_" in atom.ffType:
			aNo_epoxide[0].append(atom.aNo)

	# the head epoxide without additional O: will be displayed as group 2
	aNo_epoxide.append(aNo_epoxide[0][:])	# epoxide head without additional O
	aNo_epoxide.append([aNo_start])	# additional O in head group
	aNo_epoxide[2].remove(aNo_start)

	# the C in tail O
	for atom in bgf1.a:
		d = len(bgftools.getShortestPath(bgf1, aNo_epoxide[1][0], atom.aNo))
		if d == 3 and "C_" in atom.ffType:
			aNo_epoxide[1].append(atom.aNo)
			aNo_end = atom.aNo	# ending point for ethylene group

	# add connected H for each epoxide group
	for i, l in enumerate(aNo_epoxide):
		for ano in l:
			atom1 = bgf1.getAtom(ano)
			# looping over all atoms in bgf
			for atom2 in bgf1.a:
				if bgf.is_bonded(atom1, atom2) and atom2.is_hydrogen():
					aNo_epoxide[i].append(atom2.aNo)
		aNo_epoxide[i].sort()
		for j in aNo_epoxide[i]:
			aNo_identified.append(j)

	#!print(aNo_epoxide)
	#!print(aNo_identified)
	#!print("")


	### identify ethyl groups: something is weird in the algorithm but looks working well 20141031
	pointer = aNo_start;
	aNo_ethylene = []; stack = [];
	pointer_atom = bgf1.getAtom(pointer)
	for i in pointer_atom.CONECT:
		if not i in aNo_identified:
			stack.append(i)
		
	while 1:
		if pointer == aNo_end:
			break;

		pointer_atom = bgf1.getAtom(pointer)
		for i in pointer_atom.CONECT:
			if not i in aNo_identified:
				stack.append(i)
		
		#!print("stack: " + str(stack))
		#!print("aNo_ethylene: " + str(aNo_ethylene))
		#!print("aNo_ethylenes: " + str(aNo_ethylenes))

		try:
			pointer = stack.pop()
		except IndexError:
			break;

		if "O_" in pointer_atom.ffType and len(aNo_ethylene) > 1:
			aNo_ethylenes.append(aNo_ethylene)
			aNo_ethylene = [];

		aNo_identified.append(pointer)
		aNo_ethylene.append(pointer)
	#!print(aNo_ethylenes)
	aNo_groups = aNo_epoxide + aNo_ethylenes


	### add charges for each group
	# epoxide group
	chgFormat = "{0:2.5f}"
	print("\n\n***** Charges *****")
	print("** Epoxides **")
	for i, l in enumerate(aNo_epoxide):
		charge = 0.0;
		for ano in l:
			atom = bgf1.getAtom(ano)
			charge += atom.charge

		if i == 0:
			print("Epoxide head group charge: " + chgFormat.format(charge))
		elif i == 1:
			print("Epoxide tail group charge: " + chgFormat.format(charge))
		elif i == 2:
			print("Epoxide head group without additional O charge: " + chgFormat.format(charge))
		elif i == 3:
			print("Additional O charge: " + chgFormat.format(charge))
		else:
			nu.die("Something wrong with printing epoxide group charge.")
		
	print("\n")

	# ethylene group
	print("** Ethylenes **")
	charge_H_all = []; charge_C_all = []; charge_O_all = [];
	for i, l in enumerate(aNo_ethylenes):
		charge = 0.0;
		charge_H = []; charge_C = []; charge_O = [];
		for ano in l:
			atom = bgf1.getAtom(ano)
			charge += atom.charge
			if "C_" in atom.ffType:
				charge_C.append(atom.charge)
				charge_C_all.append(atom.charge)
			elif "H_" in atom.ffType:
				charge_H.append(atom.charge)
				charge_H_all.append(atom.charge)
			elif "O_" in atom.ffType:
				charge_O.append(atom.charge)
				charge_O_all.append(atom.charge)


		meanC, stdC = nu.meanstdv(charge_C)
		meanH, stdH = nu.meanstdv(charge_H)
		meanO, stdO = nu.meanstdv(charge_O)
		print("Ethylene group " + str(i) + " charge: " + chgFormat.format(charge))
		print(" -- C: mean " + chgFormat.format(meanC) + ", stdev " + chgFormat.format(stdC))
		print(" -- H: mean " + chgFormat.format(meanH) + ", stdev " + chgFormat.format(stdH))
		print(" -- O: mean " + chgFormat.format(meanO))
		print("")
		
	meanC, stdC = nu.meanstdv(charge_C_all)
	meanH, stdH = nu.meanstdv(charge_H_all)
	meanO, stdO = nu.meanstdv(charge_O_all)
	print("Averaging for all ethylene groups")
	print(" -- C: mean " + chgFormat.format(meanC) + ", stdev " + chgFormat.format(stdC))
	print(" -- H: mean " + chgFormat.format(meanH) + ", stdev " + chgFormat.format(stdH))
	print(" -- O: mean " + chgFormat.format(meanO) + ", stdev " + chgFormat.format(stdO))
	print("")

	# total charge
	charge = 0;
	for atom in bgf1.a:
		charge += atom.charge

	print("Total charge: " + chgFormat.format(charge))


	### save output to bgf
	bgf1.saveBGF(bgffile[:-4] + ".out.bgf")
	

	return 1

main()

