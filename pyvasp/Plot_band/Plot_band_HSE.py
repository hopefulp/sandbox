import sys,time
import matplotlib.pyplot as plt
import numpy as np
from pylab import *

version = 20130115

nt = time.localtime()
now_time = "%s_%s_%s_%s%s%s" % (nt[0],nt[1],nt[2],nt[3],nt[4],nt[5])
usage = "Usage: %s [systemlabel.band] [E min] [E max]" % sys.argv[0]
foottext = '\n Thank you\n## Kim,Hyo Seok (KAIST) <softmax1986@kaist.ac.kr>'

print "## Plotting band structure..."
print "## Version : %s \n" % version

if len(sys.argv) != 4:
    print usage
    print foottext
    sys.exit(1)

if len(sys.argv) == 4:
    filename = str(sys.argv[1])
    ymn = float(sys.argv[2])
    ymx = float(sys.argv[3])

f = open(filename, "r")
E_fermi = []
Data_of_bands = []

list_lines = []
for line in f.readlines():
	list_line = line.split()
	if list_line == []:
		continue
	else:	
		list_lines.append(list_line)
E_fermi.append(float(list_lines[0][0]))

new_list_lines = []
for x in range(1,len(list_lines)):
	temp_list = []
	for y in range(0,len(list_lines[x])):
		temp = float(list_lines[x][y])
		temp_list.append(temp)
	new_list_lines.append(temp_list)

temp_list = []
for x in range(0,len(new_list_lines)):
	temp_list.append(new_list_lines[x][0])
num_band_line = int(max(temp_list))

div = int(len(new_list_lines)/num_band_line)

count = 0
bands_line = []
for x in range(0,div):
	bands = []
	for y in range(1+count,num_band_line+1+count):
		bands.append(new_list_lines[y][1])
	count = count + num_band_line +1
	bands_line.append(bands)
	if count >= len(new_list_lines) :
		break

for x in range(0,len(bands_line[0])):
	Data_of_bands.append([])
#print Data_of_bands
for x in range(0,len(bands_line)):
	for y in range(0,len(bands_line[x])):
		#print bands_line[x][y]
		Data_of_bands[y].append(bands_line[x][y])
#print Data_of_bands[0]
		
x_val= []
for x in range(1,len(Data_of_bands[0])+1):
	x_val.append(x)

fig = plt.figure()
fig1 = fig.add_subplot(111)
#x_ticks([])

ymin = ymn
ymax = ymx
xmax = max(x_val)
xmin = min(x_val)

for x in range(0,len(Data_of_bands)):
	list_x = x_val
	list_y = Data_of_bands[x]
	array_x = np.array(list_x); array_y = np.array(list_y)
	fig1.plot(array_x, array_y, color = 'black', linewidth = '2')
	del list_y; del list_x

fig1.set_xlim(xmin, xmax)
fig1.set_ylim(ymin,ymax)
plt.show()
