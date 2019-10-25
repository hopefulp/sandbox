#!/home/noische/program/python27/bin/python
"""
PEI_NNstrain.py
Original: Aug 23 2011 In Kim
"""

# Python Modules
import sys
import os
import string
import getopt
import time

# Custom Modules
import numpy
import bgf
import bgftools

# Globals
version = '110823'


def pei_nnstrain(bgf_file, silent=True):
	"""
	"""
	# Initialization
	output = "";
	c1 = ""; c5 = ""; n = ""; m = ""; temp = ""; dist = 0; l_dist = []; l_intra_dist = []; l_inter_dist = []; nnpath = ""; count = 0;

	# Open BGF
	if isinstance(bgf_file, bgf.BgfFile):
		myBGF = bgf_file
	else:
		if not silent: print("opening bgf file.. " + str(bgf_file))
		myBGF = bgf.BgfFile(bgf_file)

	# branch-off bgf
	#myBGF2 = bgftools.getBackbone(myBGF)

	output = "n.aNo" + "\t" + "c1.aNo" + "\t" + "c5.aNo" + "\t" + "m.aNo" + "\t" + "dist" + "\t" + "nnpath"
	print(output)

	inter = 0; intra = 0;
	pri_pri = 0; pri_sec = 0; sec_sec = 0; unknown_xlnk_type = 0; saturated_pri = 0;
	n_O_atom = 0;

	for atom in myBGF.a:
		# count the number of oxygen atoms
		if "O" in atom.aName:
			n_O_atom += 1;

		# find XLK C1
		# find XLK C5 connected with C1
		### same rNo, same rName, usually C1.aNo + 3
		if "XLK" in atom.rName and "C1" in atom.aName:
			c1 = atom
			c5 = myBGF.getAtom(myBGF.a2i[c1.aNo + 3])
			count += 1;
		else:
			continue;

		# find N connected to C1: n
		for number in c1.CONECT:
			temp = myBGF.getAtom(number)
			if "N" in temp.aName:
				if bgf.is_bonded(temp, c1):
					n = temp

		# determine what kind of amine n was
		d_amine_info = dict()	# contains what was the original amine group
		connected_ano = []; n_carbon = 0; n_x_carbon = 0;
		connected_ano = n.CONECT
		for aNo2 in connected_ano:
			testatom = myBGF.getAtom(aNo2)
			if "C_" in testatom.ffType:
				if "X" not in testatom.rName:
					n_carbon += 1
				else:
					n_x_carbon += 1
		if n_carbon == 1 and n_x_carbon == 1:
			d_amine_info[n.aNo] = "PRI"
		elif n_carbon == 2 and n_x_carbon == 1:
			d_amine_info[n.aNo] = "SEC"
		elif n_carbon == 1 and n_x_carbon == 2:
			d_amine_info[n.aNo] = "PRS"	# saturated 1" amine
			saturated_pri += 1;
		else:
			d_amine_info[n.aNo] = "???"

		# find N connected to C5: m
		for number in c5.CONECT:
			temp = myBGF.getAtom(number)
			if "N" in temp.aName:
				if bgf.is_bonded(temp, c5):
					m = temp

		# determine what kind of amine m was
		connected_ano = []; n_carbon = 0; n_x_carbon = 0;
		connected_ano = m.CONECT
		for aNo2 in connected_ano:
			testatom = myBGF.getAtom(aNo2)
			if "C_" in testatom.ffType:
				if "X" not in testatom.rName:
					n_carbon += 1
				else:
					n_x_carbon += 1
		if n_carbon == 1 and n_x_carbon == 1:
			d_amine_info[m.aNo] = "PRI"
		elif n_carbon == 2 and n_x_carbon == 1:
			d_amine_info[m.aNo] = "SEC"
		elif n_carbon == 1 and n_x_carbon == 2:
			d_amine_info[m.aNo] = "PRS"	# saturated 1" amine
			saturated_pri += 1;
		else:
			d_amine_info[m.aNo] = "???"

		## check the cross-linking type: pri_pri, pri_sec, sec_sec
		if "PR" in d_amine_info[n.aNo] and "PR" in d_amine_info[m.aNo]:
			pri_pri += 1;
		elif "PR" in d_amine_info[n.aNo] and "SE" in d_amine_info[m.aNo]:
			pri_sec += 1;
		elif "SE" in d_amine_info[n.aNo] and "PR" in d_amine_info[m.aNo]:
			pri_sec += 1;
		elif "SE" in d_amine_info[n.aNo] and "SE" in d_amine_info[m.aNo]:
			sec_sec += 1;
		else:
			unknown_xlnk_type += 1;

		# check that they are connected
		#nnpath = bgftools.getShortestPath(myBGF, n.aNo, m.aNo)
		
		# distance
		dist = bgf.distance(n, m)
		l_dist.append(dist)

		# intra/inter crosslinking check
		dist = 0; output = ""
		if n.rNo != m.rNo:
			inter += 1
			dist = bgf.distance(n, m)
			l_inter_dist.append(dist)
			output += "inter\t"
		else:
			intra += 1
			dist = bgf.distance(n, m)
			l_intra_dist.append(dist)
			output += "intra\t"

		# output: N# N# distance len(nnpath)
		output += str(n.aNo) + "\t" + str(c1.aNo) + "\t" + str(c5.aNo) + "\t" + str(m.aNo) + "\t" + str(d_amine_info[n.aNo]) + "\t" + str(d_amine_info[m.aNo])  + "\t" + "{0:8.5f}".format(dist) + "\t" + str(n.rNo) + "\t" + str(m.rNo)
		print(output)

	# output: overall
	print("\nNumber of cross-linkers: " + str(count))
	print("Number of Oxygen atoms: " + str(n_O_atom))
	if count != n_O_atom:
		print("WARNING: the number of cross-linker and O atom number is different!")
	print("Number of PRI-PRI cross-linking: " + str(pri_pri))
	print("Number of PRI-SEC cross-linking: " + str(pri_sec))
	print("Number of SEC-SEC cross-linking: " + str(sec_sec))
	print("Number of saturated PRI: " + str(saturated_pri))
	print("")

	### overall cross-linking
	# output: average
	print("\nAverage of overall N-N distances with cross-linker: " + "{0:8.5f}".format(numpy.average(l_dist)))

	# output: historgram
	print("\nHistogram of overall N-N distances with cross-linker (histogram):")
	bins = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0]
	a_hist, temp = numpy.histogram(l_dist, bins)
	for index, item in enumerate(a_hist):
		print(str(temp[index]) + "\t" + str("#"*item))

	print("\nHistogram of overall N-N distances with cross-linker (numbers):")
	bins = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0]
	a_hist, temp = numpy.histogram(l_dist, bins)
	for index, item in enumerate(a_hist):
		print(str(temp[index]) + "\t" + str(item))


	### intra cross-linking
	# output: average
	print("\nAverage of intra N-N distances with cross-linker: " + "{0:8.5f}".format(numpy.average(l_intra_dist)))

	# output: historgram
	print("\nHistogram of intra N-N distances with cross-linker (histogram):")
	bins = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0]
	a_hist, temp = numpy.histogram(l_intra_dist, bins)
	for index, item in enumerate(a_hist):
		print(str(temp[index]) + "\t" + str("#"*item))

	print("\nHistogram of intra N-N distances with cross-linker (numbers):")
	bins = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0]
	a_hist, temp = numpy.histogram(l_intra_dist, bins)
	for index, item in enumerate(a_hist):
		print(str(temp[index]) + "\t" + str(item))


	### inter cross-linking
	# output: average
	print("\nAverage of inter N-N distances with cross-linker: " + "{0:8.5f}".format(numpy.average(l_inter_dist)))

	# output: historgram
	print("\nHistogram of inter N-N distances with cross-linker (histogram):")
	bins = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0]
	a_hist, temp = numpy.histogram(l_inter_dist, bins)
	for index, item in enumerate(a_hist):
		print(str(temp[index]) + "\t" + str("#"*item))

	print("\nHistogram of inter N-N distances with cross-linker (numbers):")
	bins = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0]
	a_hist, temp = numpy.histogram(l_inter_dist, bins)
	for index, item in enumerate(a_hist):
		print(str(temp[index]) + "\t" + str(item))

	print("\nNumber of intra crosslinkings: " + str(intra))
	print("\nNumber of inter crosslinkings: " + str(inter))

	##### End of pei_nnstrain


if __name__ == '__main__':

	option = ""; args = ""; bgf_file = ""; ff_file = ""; probability = 0; out_file = "";
	usage = """
Usage: 	
"""

	if len(sys.argv) < 2:
		print(usage)
		sys.exit(0)

	# Defaults
	safemode = False; silent = True;

	options, args = getopt.getopt(sys.argv[1:], 'hb:', ['help','bgf='])
	for option, value in options:
		if option in ('-h', '--help'):
			print usage; sys.exit(0)
		elif option in ('-b', '--bgf'):
			bgf_file = value
		elif option in (''):
			print usage; sys.exit(0)

	pei_nnstrain(bgf_file, silent=False)
