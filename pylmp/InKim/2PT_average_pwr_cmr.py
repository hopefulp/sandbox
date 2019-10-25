#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import pickle
import sys
import numpy as np

f = open('pwr_spectrum.pickle')
result = pickle.load(f)
t = result.keys()
t.sort()
last_t = t[-100:]
x = np.array(result[last_t[-1]]['freq(cm-1)'])

sum = np.zeros(x.shape)
for i in last_t:
    sum += result[i]['PWRcmr'][2]
sum /= len(last_t)

for index, i in enumerate(x):
    print "%8.5f %8.5f" % (x[index], sum[index])
