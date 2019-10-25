#!/opt/applic/epd/bin/python
import cPickle as pkl
import sys
import numpy as np

"""
when I get the dipole distribution data by ACETONE_getDipoleDistr.py, water-acetone-layer_final2-0_100K_10ns.dipoleDistr.pickle will be generated.
This file is to average for timesteps contained in that pickle file.

notice that this is only for NVT cases

pickle file description:

timestep - residue - [ [bin angle magnitude] , ... ]

"""

f = open(sys.argv[1])
result = pkl.load(f)

timesteps = result.keys()
timesteps.sort()

residues = result[timesteps[0]].keys()
residues.sort()

bins = [];
for i in result[timesteps[0]][residues[0]]:
	bins.append(i[0])	# nvt

# allocate for the numpy array

temp = [];
for r in residues:
	temp2 = np.zeros((len(bins), 3))
	for t in timesteps:
		temp2 += result[t][r]
	temp2 /= len(timesteps)
	temp.append([r, temp2])

# print
output = "";
for i in temp:
	print(str(i[0]))
	for j in i[1]:
		print(str(j[0]) + '\t' + str(j[1]) + '\t' + str(j[2]))
