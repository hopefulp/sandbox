#!/home/joonho/anaconda3/bin/python

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import MaxNLocator

def file_len(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i+1

name = 'neb.dat'
Tlines = file_len('%s' % name)

data_all = []
f = open('%s' % name, 'r')
for i in range(Tlines):
    line = f.readline()
    words = line.split()
    data_all.append(words)
f.close()

image  = []; dist  = []
energy = []; force = []

# Data gathering
for i in range(len(data_all)):
    image.append(int(data_all[i][0]))
    dist.append(float(data_all[i][1]))
    energy.append(float(data_all[i][2]))
    force.append(float(data_all[i][3]))

maxindex = energy.index(max(energy))
minindex = energy.index(min(energy))

# Visualization
fig, ax1 = plt.subplots(figsize=(8,7.5))
ax2 = ax1.twinx()

lines1 = ax1.plot(image, energy, color='black', label='Energy')
lines2 = ax2.plot(image, force, color='blue', label='Force')

ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
ax1.set_xlabel('Reaction Coordinate', fontsize=15)
ax1.set_ylabel('Energy [eV]', fontsize=15)
ax2.set_ylabel('Force [eV / $\AA$]', fontsize=15, color='blue')

lines = lines1+lines2
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc=0, fontsize=12)
ax1.text(image[maxindex]+0.2, 0.6*energy[maxindex], "Barrier=%4.2f eV" % energy[maxindex], fontsize=15)
ax1.vlines(image[maxindex],energy[minindex],energy[maxindex],color='black',linestyles='--')
ax2.hlines(0,image[0], image[-1],color='blue',linestyles='--')
plt.tight_layout()
plt.savefig('neb_result.png', dpi=300, bbox_inches ="tight")
plt.show()


