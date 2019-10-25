#!/usr/bin/python2.7
#heejin
import sys
import os
import re
import math
import subprocess
import operator
import itertools

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

#
if not os.path.isfile('EIGENVAL'):
	print "EIGENVAL is not found."
	sys.exit()
if not os.path.isfile('OUTCAR'):
	print "OUTCAR is not found."
	sys.exit()

if len(sys.argv) is 2:
	eshift = float(sys.argv[1])
else:
	eshift = 0

#read eigenval
eigenfile = open('EIGENVAL')

line = eigenfile.readline()
line = eigenfile.readline()
line = eigenfile.readline()
line = eigenfile.readline()
line = eigenfile.readline()
line = eigenfile.readline()
nband = line.split()[2]
npoint = line.split()[1]
eigenfile.close()

ef = float(getsh("grep E-fermi OUTCAR")[0].strip().split()[2])

#print "%10.6f  %10s  %10s" % (ef, nband, npoint)

totset = []
for i in range (int(nband)+1): #i:0-6
	totset.append([0]*(int(npoint)+1)) #col:41

for i in range(1,int(nband)+1): #i:1-6
	grpband = "grep ' %d     ' EIGENVAL" % i
	res = getsh(grpband)
	for j in range(1,int(npoint)+1): #j:1-40
		#totset[0][0] = 0
		totset[0][j] = str(j) + ' ' + str(j)
		spres = res[j-1].strip().split() #res[0]-res[39]
		outline = "%10.6f %10.6f" % (float(spres[1])-ef+eshift, float(spres[2])-ef+eshift)
		totset[i][j] = outline #i:1-6
#print totset

revset = []
for i in range (int(npoint)+1): #i:0-40
	revset.append([0]*(int(nband)+1)) #col:7

revset[0][0] = ' '
for i in range(0,int(nband)+1): #i:0-6
	for j in range(1,int(npoint)+1): #j:1-40
		revset[0][i] = 'k' + str(i) + 'up' + ' ' + 'k' + str(i) + 'down'
		revset[j][i] = str(totset[i][j]) #revset[1][0]-revset[40][6]

sumset = []
for i in revset:
	temp = ''
	for j in i:
		temp += str(j) + ' '
	sumset.append(temp)

for i in range(int(npoint)+1):
	print sumset[i]
