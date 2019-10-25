#!/opt/applic/epd/bin/python
import cPickle as pkl
import numpy as np

"""
parse pickle which contains dipole moment and get acetone dipole distribution at interfaces (defined in x)
"""
f = open('water-acetone-layer_final2-0_100K_10ns_Eq50ns.dipole.histo.pickle')
a = pkl.load(f)
bins = np.arange(-1, 1, 0.01)

x = [39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0]
y = [a['ACT'][i] for i in x]
n = [[] for i in y]

for index, i in enumerate(y):
	n[index], bin_edges = np.histogram(i, bins)

output = "";
for index1, i in enumerate(n):
	for j in i:
		output += str(j) + '\t'
	output += '\n'

print(output)

