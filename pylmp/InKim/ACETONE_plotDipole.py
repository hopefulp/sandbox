#!/opt/applic/epd/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

import sys
import os
import math
import pickle

### initialize
fig = plt.figure()
ax = fig.add_subplot(111)
axcolor = 'lightgoldenrodyellow'

### data read
myRESULT = open(sys.argv[1])
result = pickle.load(myRESULT)
print("Load successful")


### REMARK: there's a problem that n and r cannot be loaded in the eventhandler.
###         so add a global value and set when the r changes.
# keywords
def get_k():
	global curr_k
	return curr_k
def set_k(val):
	global curr_k
	curr_k = val

# timesteps
def get_t():
	global curr_t
	return curr_t
def set_t(val):
	global curr_t
	curr_t = val



### get a data
# result: timestep - MU, ABSMU, ANGELS, NATOMS
timesteps = result.keys()
timesteps.sort()
t = timesteps[0]	# the first timestep
n_t = len(timesteps) - 1	# number of timesteps

keywords = result[t].keys()
k = keywords[0]
set_k(k)

residues = result[t][k].keys()

data = [];
for i in residues:
	data.append(result[t][k][i])

x = result[t][k]
### draw histogram
bins = np.arange(0, 181, 1);	# 0 ~ 2pi
#n, bin_edges = np.histogram(data, bins, density=True)
n, bin_edges = np.histogram(data, bins)

count = 0;
for i in n:
	count += i;
print(count)
#### normalize
#m = [];
#for index, i in enumerate(n):
#	if index != 0:
#		m.append(n[index] / np.sin(np.deg2rad(bins[index])))
#	else:
#		m.append(n[index] / np.sin(np.deg2rad(0.0001)))

# truncate the last 180
bins = bins[:-1]

#l, = plt.plot(bins, m)
l, = plt.plot(bins, n)

ax_timestep = plt.axes([0.25, 0.03, 0.65, 0.03], axisbg=axcolor)
s_timestep = Slider(ax_timestep, 'Timestep', 0, n_t, valinit=t, valfmt='%3d')

def update(val):
	data = [];
	for i in residues:
		data.append(result[timesteps[int(val)]][k][i])

	bins = np.arange(0, 181, 1);	# 0 ~ 2pi
	#n, bin_edges = np.histogram(data, bins, density=True)
	n, bin_edges = np.histogram(data, bins)
	bins = bins[:-1]

	count = 0;
	for i in n:
		count += i;
	print(count)

	#m = [];
	#for index, i in enumerate(n):
	#	if index != 0:
	#		m.append(n[index] / np.sin(np.deg2rad(bins[index])))
	#	else:
	#		m.append(n[index] / np.sin(np.deg2rad(0.0001)))

	#l.set_ydata(m)
	l.set_ydata(n)
	plt.draw()

s_timestep.on_changed(update)

plt.show()
