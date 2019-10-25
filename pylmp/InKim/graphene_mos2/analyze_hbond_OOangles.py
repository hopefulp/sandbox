#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys
import numpy as np
import cPickle as pkl

pkl_file = sys.argv[1]

with open(pkl_file, 'rb') as f:
    angles = pkl.load(f)

t = sorted(angles.keys())

Angles = []
for a in angles:
    for i in angles[a]:
        Angles.append(i)

delta_angle = 1.0
bin_angle = np.arange(0, 180 + delta_angle, delta_angle)

distr_angle, _ = np.histogram(Angles, bins=bin_angle)
norm_distr_angle = distr_angle/float(distr_angle.sum())

output = ''
for i, j in zip(_[1:], norm_distr_angle):
    output += "%5.3f" % (i - delta_angle/2.0) + "%8.5f" % j + "\n"

print output
