#!/home/noische/python

from VASP import *
import numpy as np
from numpy import sqrt
import bgftools

m = VASP()
#h = m._latticeVectors = [[5.0, 0, 0], [5*sqrt(2)/2, 5*sqrt(2)/2, 0], [0, 0, 5.0]]
h = m._latticeVectors = [[5.0, 0, 0], [5, 5, 0], [0, 0, 5.0]]
hinv = m._inv_latticeVectors = np.linalg.inv(m._latticeVectors)
#dim = [1.0, 1.0, 2.0]
#points = [[0, 0, 0], [1, 1, 2], [1, 0, 2], [0.5, 0.5, 1.5]]
points = [[2, 1, 0], [2, 3, 0]]
for i in points:
	for j in points:
		if i==j: continue;
		print("vasp dist: " + str(i) + str(j) + str(m.distance(i, j, h, hinv)))
		#print("pbc_dist: " + str(i) + str(j) + str(bgftools.pbc_dist(i, j, dim)))
		print("")

print(sqrt(5))
print(sqrt(10))
