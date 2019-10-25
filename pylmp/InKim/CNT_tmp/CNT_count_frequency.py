#!/home/noische/python

"""
A script for the project /qcfs/noische/research/CNT-fluctuation

Count the frequency from the water count profile, i.e., 
fr  n_wat   r   h
1       3   4.022     8.8
2       3   4.172     8.8
3       3   4.089     8.9
"""

import sys
import numpy as np
import nutils as nu

profile = sys.argv[1]
f = open(profile)

f.readline()    # pass the title 

n_wat = []
while 1:
    try:
        line = f.readline()

        if not line:
            break;

        line = line.replace("\n", "")
        line = line.split()
        n_wat.append(int(line[1]))

    except KeyboardInterrupt:
        nu.die("Keyboard Break.. Exiting.")

f.close()

n_wat = np.array(n_wat)
unique, counts = np.unique(n_wat, return_counts=True)

#print np.asarray((unique, counts)).T
for index, i in enumerate(unique):
    print "%d\t%d" % (unique[index], counts[index])

