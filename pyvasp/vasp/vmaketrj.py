#!/usr/bin/python2.7
#heejin
import sys
import os
import re
import math
import subprocess
import operator

#usage description: nstep
if len(sys.argv) < 3:
	print "Usage: [step] [atom name] [outname] [last image:optional]"
	sys.exit()

nstep=sys.argv[1]
targetname=sys.argv[2]
outhead=sys.argv[3]

if len(sys.argv) == 5:
	laststep = sys.argv[4]

if not os.path.isfile('POSCAR'):
	print "POSCAR is not found."
	sys.exit()
if not os.path.isfile('XDATCAR'):
	print "XDATCAR is not found."
	sys.exit()

#copy XDATCAR
os.system('cp XDATCAR XDATCAR-trj')
os.system('sed -i 3d XDATCAR-trj')
os.system('sed -i 3d XDATCAR-trj')
os.system('sed -i 3d XDATCAR-trj')

#read poscar: species, nspecies
poscarfile = open('POSCAR')

species = []
line = poscarfile.readline()
species = line.split()

for i in range (4):
	poscarfile.readline()

line = poscarfile.readline()
if line.split()[0].isalpha():
	line = poscarfile.readline()
	nspecies = line.split()
else:
	nspecies = line.split()
poscarfile.close()

#calc number of atoms: ntotal, ntarget
ntotal = 0
ntarget = 0
for i in range(len(nspecies)):
	if species[i] in targetname:
		ntarget = int(nspecies[i])
	ntotal = ntotal + int(nspecies[i])
if ntarget is 0:
	print targetname +" is not in POSCAR"
	sys.exit()

print "Trajectory is generating for " + str(ntarget) + " " + targetname + " atoms"

#count XDATCAR
nlxdatcar = sum(1 for line in open('XDATCAR'))
numset = (nlxdatcar - 7) / (ntotal + 1)

if (len(sys.argv) == 5 and int(laststep) < int(numset)):
	numset = int(laststep)

numtarget = str((numset / int(nstep) + 2) * ntarget)
natomstr = "  ".join(nspecies[0:-1]) + "  " + numtarget + "\n"

outputname = outhead+"_i"+str(numset)+"s"+nstep+"_POSCAR"
numframe = str(numset / int(nstep))

print " - Total " + str(numset) + " images"
print " - Writing " + numframe + " images for each " + str(nstep) + " step"
print " - Output file name is " + outputname

#write new POSCAR
poscarfile = open('POSCAR')
newfile = open(outputname,'w')
while 1:
	line = poscarfile.readline()
	if not line:
		break
	if line.strip() == '':
		break
	if ' ' + str(ntarget) in line:
		line = poscarfile.readline()
		newfile.write(natomstr)
	newfile.write(line)
poscarfile.close()
newfile.close()

#build script
scriptfile = open('temp.sh','w')
do0 = '#!/bin/tcsh\n'
do1 = 'set numlist = `seq 1 ' + str(nstep) + ' ' + str(numset) + '`\n'
do2 = 'foreach nl(${numlist})\n'
do3 = 'xdat2postrj.pl 1 $nl\n'
do4 = 'end\n'
do5 = 'foreach sect(`ls POSCAR*.out`)\n'
do6 = 'tail -n ' + str(ntarget) + ' $sect >> ' + str(outputname) + '\n'
do7 = 'end\n'
do8 = 'rm POSCAR*.out\n'

scriptfile.write(do0)
scriptfile.write(do1)
scriptfile.write(do2)
scriptfile.write(do3)
scriptfile.write(do4)
scriptfile.write(do5)
scriptfile.write(do6)
scriptfile.write(do7)
scriptfile.write(do8)
scriptfile.close()

os.system('chmod +x temp.sh')
os.system('temp.sh')
os.system('rm temp.sh XDATCAR-trj')
