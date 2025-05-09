#!/home/joonho/anaconda3/bin/python

from NanoCore import *
from matplotlib import cm
import time, sys,os,glob
import argparse
import siesta as s2
#label = sys.argv[1]

### siesta vs vasp

sim = 'vasp'

if sim == 'siesta':
    at = io.read_xyz('CdS2_H_term.xyz') # % label)
elif sim == 'vasp':
    at = io.read_poscar('POSCAR')
atoms = at._atoms
#print atoms
# slice atoms by z coordinates
z_coords = []; indice = []
for atom in atoms:
    if not atom[2] in z_coords: z_coords.append(atom[2])

for z in z_coords:
    temp = []
    for atom in atoms:
        if abs(z-atom[2]) < 0.01: temp.append(atom.get_serial())
    indice.append(temp)
#print indice
simobj = 0.0
### find fermilevel for siesta and vasp
if sim == 'siesta':
    a=glob.glob('*.EIG')
    a = str(a[0])
    f=open(a)
    list_lines=[]
    for line in f.readlines():
        list_lines.append(line)
    Fermi = list_lines[0]
    Fermi = Fermi.split()
    Fermi = float(Fermi[0])
elif sim == 'vasp':
    with open('DOSCAR', 'r') as f:
        for i, line in enumerate(f):
            if i == 5:
                ele = line.strip().split()
                Fermi = float(ele[3])
                break
print (Fermi)
# get pdos
Z = []; E = []
for ind in indice:
    E1, dos11, dos12 = s2.get_pdos(simobj, -15, 5, by_atom=1, atom_index=ind, broad= 0.02, npoints = 1500, label = 'CdS2_H_term')
    E = np.array(E1)
    Z.append(np.array(dos11))

absZ = np.abs(Z).T
print (len(E))
# convert to log10 values + minimum correction to avoid -INF
Z = np.log10(absZ + 10**-5)

# generate meshgrid
X, Y = np.meshgrid(np.array(z_coords), E-Fermi)

# customized figure 
import matplotlib.pyplot as plt
fig1 = plt.figure(figsize=(10,12))
plt.yticks(np.arange(-4,4,1))
levels = np.linspace(1.01*Z.min(), 0.99*Z.max(), 100)
cmap=plt.cm.get_cmap("jet")
import pylab as plb
cset = plb.contourf(X,Y,Z, levels, cmap=cmap)
plb.colorbar(cset,ticks=[-4,-3,-2,-1,0,1])
#plt.show()
fig1.savefig('pldos.png')



