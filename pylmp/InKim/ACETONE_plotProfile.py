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
curr_r = "";
curr_t = "";

### data read
myRESULT = open(sys.argv[1])
result = pickle.load(myRESULT)
print("Load successful")


### REMARK: there's a problem that n and r cannot be loaded in the eventhandler.
###         so add a global value and set when the r changes.
def set_r(val):
	global curr_r
	curr_r = val
def set_t(val):
	global curr_t
	curr_t = val

def get_r():
	global curr_r
	return curr_r
def get_t():
	global curr_t
	return curr_t



### get a data
# result: timestep - residue - HIST,BIN_EDGES
timesteps = result.keys()
timesteps.sort()
t = timesteps[0]	# the first timestep
n_t = len(timesteps) - 1	# number of timesteps

residues = result[t].keys();
r = residues[0]	# ANGLES
set_r(r)

x = result[t][r]['BIN_EDGES'][:-1]
y = result[t][r]['HIST']

l, = plt.plot(x, y)


# timestep slider
ax_timestep = plt.axes([0.25, 0.03, 0.65, 0.03], axisbg=axcolor)
s_timestep = Slider(ax_timestep, 'Timestep', 0, n_t, valinit=t, valfmt='%5d')

def t_update(val):
	n = timesteps[int(math.ceil(val))]

	set_t(n)
	r = get_r();

	y = result[n][r]['HIST']
	l.set_ydata(y)
	plt.draw()
s_timestep.on_changed(t_update)


# residue selector
rax = plt.axes([0.025, 0.5, 0.15, 0.15], axisbg=axcolor)
radio = RadioButtons(rax, residues, active=0)

def r_update(val):
	#n = timesteps[int(math.ceil(s_timestep.val))]
	n = get_t();

	r = str(val)
	set_r(val)

	y = result[n][r]['HIST']
	l.set_ydata(y)
	plt.draw()
radio.on_clicked(r_update)


# print button
ax_print = plt.axes([0.8, 0.1, 0.1, 0.04])
button = Button(ax_print, 'Print', color=axcolor, hovercolor='0.975')
def printstep(event):
	n = get_t();
	r = get_r();

	x = result[n][r]['BIN_EDGES'][:-1]
	y = result[n][r]['HIST']

	print('TIMESTEP ' + str(n) + '\t' + str(r))
	for index, i in enumerate(x):
		print(str(i) + '\t' + str(y[index]))
button.on_clicked(printstep)
	
plt.show()
