#!/usr/bin/python2.7
#heejin
import sys
import os
import re
import subprocess
import operator
import math

def getsh(command):
	pp = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
	stdout = pp.stdout
	result = []
	while 1:
		line = stdout.readline()
		if not line:
			break
		result.append(line)
	return result

#read poscar
if os.path.isfile('CONTCAR'):
	posfile = open('CONTCAR')
elif os.path.isfile('POSCAR'):
	posfile = open('POSCAR')
else:
	print "POSCAR is not found."
	sys.exit()

#ignore head
line = posfile.readline() #species
line = posfile.readline() #scale
line = posfile.readline() #a
line = posfile.readline() #b
line = posfile.readline() #c
line = posfile.readline() #species

#count total number of ions
line = posfile.readline() #no.species
nspecies = line.split()
ntotal = 0
for i in range(len(nspecies)):
	ntotal = ntotal + int(nspecies[i])

#main
saveforce = "grep -A %d TOTAL-FORCE OUTCAR" % (ntotal + 1)
force = ['--']
force = force + getsh(saveforce)

stp = 0
for i in range(len(force)):
	line = force[i].strip()
	if '-----------------------------------' in line:
		stp = stp + 1
		forcestep = []
		for j in xrange((ntotal+3)*(stp-1)+3, ((ntotal+3)*stp)):
			lline = force[j].strip().split()
			sca = math.sqrt(float(lline[3])*float(lline[3])+float(lline[4])*float(lline[4])+float(lline[5])*float(lline[5]))
			forcestep.append(sca)
		print "%4d  %6.3f" % (stp, max(forcestep))
