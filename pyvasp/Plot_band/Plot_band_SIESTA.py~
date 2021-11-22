#!/usr/bin/env python
#------------------------------------------------------------------------
# This script plots the E-k band diagram reading the SEISTA output file #
# systemlabel.bands. If there are bands for both spin states, band for  #
# each one are plotted                                                  #
#                                                                       #
#  Script by Javier Junquera and Pablo Aguado-Puente (2007)             #
#  Updated by Hyo Seok Kim (2013)                                       #                
#------------------------------------------------------------------------

import sys,time
import matplotlib.pyplot as plt
import numpy as np
from pylab import *


version = 20130115

nt = time.localtime()
now_time = "%s_%s_%s_%s%s%s" % (nt[0],nt[1],nt[2],nt[3],nt[4],nt[5])
usage = "Usage: %s [systemlabel] [E min] [E max] [PDOS_source_1] [PDOS_sorce_2] [PDOS_sorce_3]" % sys.argv[0]
foottext = '\n Thank you\n## Kim,Hyo Seok (KAIST) <softmax1986@kaist.ac.kr>'

print "## Plotting band structure..."
print "## Version : %s \n" % version

#Check input file

if len(sys.argv) == 4:
    filename = str(sys.argv[1]+'.bands')
    ymn = float(sys.argv[2])
    ymx = float(sys.argv[3])
    DOS_file = str(sys.argv[1]+'.DOS')

elif len(sys.argv) == 5:
    filename = str(sys.argv[1]+'.bands')
    ymn = float(sys.argv[2])
    ymx = float(sys.argv[3])
    DOS_file = str(sys.argv[1]+'.DOS')
    PDOS_file_1 = str(sys.argv[4]+'.dat')

elif len(sys.argv) == 6:
    filename = str(sys.argv[1]+'.bands')
    ymn = float(sys.argv[2])
    ymx = float(sys.argv[3])
    DOS_file = str(sys.argv[1]+'.DOS')
    PDOS_file_1 = str(sys.argv[4]+'.dat')
    PDOS_file_2 = str(sys.argv[5]+'.dat')

elif len(sys.argv) == 7:
    filename = str(sys.argv[1]+'.bands')
    ymn = float(sys.argv[2])
    ymx = float(sys.argv[3])
    DOS_file = str(sys.argv[1]+'.DOS')
    PDOS_file_1 = str(sys.argv[4]+'.dat')
    PDOS_file_2 = str(sys.argv[5]+'.dat')
    PDOS_file_3 = str(sys.argv[6]+'.dat')

else:
    print usage
    print foottext
    sys.exit(1)


#Read the input file where the band structure is stored.

f = open(filename,"r")
f1 = open(DOS_file,"r")

if len(sys.argv) == 5:
    f2 = open(PDOS_file_1,"r")

if len(sys.argv) == 6:
    f2 = open(PDOS_file_1,"r")
    f3 = open(PDOS_file_2,"r")

if len(sys.argv) == 7:
    f2 = open(PDOS_file_1,"r")
    f3 = open(PDOS_file_2,"r")
    f4 = open(PDOS_file_3,"r")

#=========================================Read data==============================================#

#Data_of_Band : [[Fermi level], [K-point Min, K-point Max], [Energy min, Energy Max],[The number of obitals, The number of spin components, The number of K-points]]
#xk_points : The list of all of the  k-points
#bands : the list of eigenvalnue energy
#Special_k_point : [['Gamma, "'K'"], ......]

Data_of_Band = []
xk_points = []
bands = []
special_k_point = []

#DOS related value

DOS_energy = []
DOS_val = []
DOS_val_spin = []

#PDOS related value
PDOS_s1_val = []
PDOS_s1_val_spin = [] 
PDOS_s2_val = []
PDOS_s2_val_spin = []
PDOS_s3_val = []
PDOS_s3_val_spin = []

#Read information about DOS

DOS_atoms = f1.readline()
n = float(DOS_atoms.split()[0])
for i in f1.readlines():
    list_lines = i.split()
    DOS_energy.append(float(list_lines[0]))
    DOS_val.append(float(list_lines[1])/n)
    DOS_val_spin.append(float(list_lines[2])/n)
    
#Read information about PDOS
if len(sys.argv) == 5:
    s1_atoms = f2.readline()
    n = float(s1_atoms.split()[0])
    for i in f2.readlines():
        list_lines = i.split()
        #print list_lines[1]
        PDOS_s1_val.append(float(list_lines[1])/n)
        PDOS_s1_val_spin.append(float(list_lines[2])/n)

if len(sys.argv) == 6:
    s1_atoms = f2.readline()
    n = float(s1_atoms.split()[0])
    for i in f2.readlines():
        list_lines = i.split()
        PDOS_s1_val.append(float(list_lines[1])/n)
        PDOS_s1_val_spin.append(float(list_lines[2])/n)
    s2_atoms = f3.readline()
    n = float(s2_atoms.split()[0])
    for i in f3.readlines():
        list_lines = i.split()
        PDOS_s2_val.append(float(list_lines[1])/n)
        PDOS_s2_val_spin.append(float(list_lines[2])/n)

if len(sys.argv) == 7:
    s1_atoms = f2.readline()
    n = float(s1_atoms.split()[0])
    for i in f2.readlines():
        list_lines = i.split()
        print list_lines[1]
        PDOS_s1_val.append(float(list_lines[1])/n)
        PDOS_s1_val_spin.append(float(list_lines[2])/n)
    s2_atoms = f3.readline()
    n = float(s2_atoms.split()[0])
    for i in f3.readlines():
        list_lines = i.split()
        PDOS_s2_val.append(float(list_lines[1])/n)
        PDOS_s2_val_spin.append(float(list_lines[2])/n)
    s3_atoms = f4.readline()
    n = float(s3_atoms.split()[0])
    for i in f4.readlines():
        list_lines = i.split()
        PDOS_s3_val.append(float(list_lines[1])/n)
        PDOS_s3_val_spin.append(float(list_lines[2])/n)
    


#Read information of part of Data of bands. 
for i in range(0,4):
    line = f.readline()
    list_line=line.split()
    Data_of_Band.append(list_line)

#norb is the number of orbital in the unit cell,
#nspin is the nummer of spin components,
#nkpoint is the number of k-point,
#nlines is the number of lines per kpoints containing the obital
norb = int(Data_of_Band[3][0]) * int(Data_of_Band[3][1])
print "Number of orbitals in unit cell           =", norb
nspin = int(Data_of_Band[3][1])
print "Number of different spin polarization     =", nspin
nkpoint = int(Data_of_Band[3][2])
print "Number of k-points to represent the band  =", nkpoint
nlines = norb / 10
nrest = norb - (nlines * 10)

#Initialize the vector hoding the bands
for initi_band in range(0, norb):
    bands.append([])


#Read the bands and xk_points
for initi_k in range(1, nkpoint + 1):
    for initi_lines in range(0, nlines):
        line = f.readline()
        t = line.split()
        if initi_lines == 0:
            xk_points.append(float(t[0]))
            for it in range(1, len(t)):
                initi_band = it - 1
                bands[initi_band].append(float(t[it])-float(Data_of_Band[0][0]))
        else:
            for it in range(0, len(t)):
                initi_band = initi_lines * 10 + it
                bands[initi_band].append(float(t[it])-float(Data_of_Band[0][0]))
    if nrest > 1:
        line = f.readline()
        t = line.split()
        for it in range(0,len(t)):
            initi_band = nlines*10+it
            bands[initi_band].append(float(t[it])-float(Data_of_Band[0][0]))


line = f.readline()
t = line.split()
nlabels = int(t[0])

#Read special k_point information
for s in range(0,int(t[0])):
    line = f.readline()
    t = line.split()
    special_k_point.append(t)


    
#====================== Plotting E-k diagram (Using Matplotlib)===========================#



#Set range of x, y axis
Fermi_Energy = float(Data_of_Band[0][0])
ymin = min(bands[0])
ymin = ymn
ymax = ymx
xmin = float(Data_of_Band[1][0])
xmax = float(Data_of_Band[1][1])

print ymin, ymax
#Drawing structure of figure
fig = plt.figure(1,figsize = (16,9))
fig1 = fig.add_subplot(131)
fig1.set_ylabel(r'$E$ $-$ $E_f$   $[eV]$', fontsize=20)
fig1.set_xlabel('')
fig1.set_xlim(xmin, xmax)
fig1.set_ylim(ymin, ymax)
xticks([])
fig2 = fig.add_subplot(132)
fig2.set_ylabel('')
fig2.set_xlabel('DOS',fontsize=20)
fig2.set_xlim(0,0.2)
fig2.set_ylim(ymin,ymax)

#Plot Special K points 
num_k = len(special_k_point)
k_point = []
name_k_point = []
for i_kp in range(0, num_k):
    k_point.append(float(special_k_point[i_kp][0]))
    name_k_point.append(special_k_point[i_kp][1])
for j_kp in range(0,num_k):
    fig1.axvline(x=k_point[j_kp], ymin=ymin, ymax=ymax, color ='r')
    fig1.text(k_point[j_kp], ymin-0.2, name_k_point[j_kp], size = 13, horizontalalignment='center',verticalalignment='center')

#Plot relative Fermilevel

fig1.axhline(0,color = 'blue',lw=1,linestyle=':')
fig2.axhline(0,color = 'blue',lw=1,linestyle=':')

#plot band structure(Spin consideration)
if nspin > 1:
    for i_plot in range(0, norb/nspin-1):
        list_X = xk_points
        list_Y = bands[i_plot]
        array_X = np.array(list_X); array_Y = np.array(list_Y)
        fig1.plot(array_X, array_Y, color = 'black', linewidth = '2')
        del list_X; del list_Y

    for k_plot in range(norb/2,norb):
        list_x = xk_points
        list_y = bands[k_plot]
        array_x = np.array(list_x); array_y = np.array(list_y)
        fig1.plot(array_x, array_y, color = 'black', linewidth = '2',linestyle ='--')
        del list_y; del list_x
else:
    for i_plot in range(0,norb):
        list_x = xk_points
        list_y = bands[i_plot]
        array_x = np.array(list_x); array_y = np.array(list_y)
        fig1.plot(array_x, array_y, color = 'black', linewidth = '2')
        del list_y; del list_x

#plotting DOS calculation
array_DOS = np.array(DOS_val)
array_DOS_spin = np.array(DOS_val_spin)
list_E = []
for i in range(0,len(DOS_energy)):
    En = DOS_energy[i] - Fermi_Energy
    list_E.append(En)
array_DE = np.array(list_E)
#print len(array_DE)
fig2.plot(array_DOS,array_DE, color = 'black', linewidth = '2')
fig2.plot(array_DOS_spin, array_DE, color = 'black', linewidth = '2',linestyle ='--')
labels = yticks([])

#ploting PDOS calculation
if len(sys.argv) > 4 :
    fig3 = fig.add_subplot(133)
    fig3.set_ylabel('')
    fig3.set_xlabel('PDOS',fontsize=20)
    fig3.set_xlim(0,0.2)
    fig3.set_ylim(ymin,ymax)
    fig3.axhline(0,color = 'blue',lw=1,linestyle=':')
    labels = yticks([])
    if len(sys.argv) >= 5 :
        array_PDOS_1 = np.array(PDOS_s1_val)
        array_PDOS_1_spin = np.array(PDOS_s1_val_spin)
        fig3.plot(array_PDOS_1,array_DE, color = 'black', linewidth = '2')
        fig3.plot(array_PDOS_1_spin,array_DE, color = 'black', linewidth = '2',linestyle ='--')
        if len(sys.argv) >= 6 :
            array_PDOS_2 = np.array(PDOS_s2_val)
            array_PDOS_2_spin = np.array(PDOS_s2_val_spin)
            fig3.plot(array_PDOS_2, array_DE, color = 'red', linewidth = '2')
            fig3.plot(array_PDOS_2_spin,array_DE, color = 'red', linewidth = '2',linestyle ='--')
            if len(sys.argv) >= 7 :
                array_PDOS_3 = np.array(PDOS_s3_val)
                array_PDOS_3_spin = np.array(PDOS_s3_val_spin)
                fig3.plot(array_PDOS_3,array_DE, color = 'yellow', linewidth = '2')
                fig3.plot(array_PDOS_3_spin,array_DE, color = 'yellow', linewidth = '2',linestyle ='--')
    
plt.show()

