#!/usr/bin/python2.7
#heejin
import sys
import os
import re
import math

def dotproduct(v1, v2):
	return sum((a*b) for a, b in zip(v1, v2))
def crossproduct(a, b):
	crs = [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]
	return crs
def length(v):
        return math.sqrt(dotproduct(v, v))
def xuniqueCombinations(items, n):
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for cc in xuniqueCombinations(items[i+1:],n-1):
                yield [items[i]]+cc\

if len(sys.argv) < 2:
	print "Usage: [file1] [file2] [cutoff]"
	sys.exit()

filename1 = sys.argv[1]
filename2 = sys.argv[2]
cutoff = float(3)

#if len(sys.argv) == 3:
#	cutoffr = "3"
#if len(sys.argv) == 4:
#	cutoffr = sys.argv[3]

if not os.path.isfile(filename1):
	print inifile + " is not found."
	sys.exit()

if not os.path.isfile(filename2):
	print finfile + " is not found."
	sys.exit()

#if cutoffr.isdigit():
#	cutoff = float(cutoffr)
#else:
#	cutoff = float(2.3)
#	print "cutoff is wrong."
#	sys.exit()


file1 = open(filename1)
file2 = open(filename2)

line1 = file1.readline() #species
line2 = file2.readline()
species = line1.split()

line1 = file1.readline() #scale
line2 = file2.readline()

line1 = file1.readline() #a
parse = line1.split()
ax1 = float(parse[0])
ay1 = float(parse[1])
az1 = float(parse[2])

line2 = file2.readline()
parse = line2.split()
ax2 = float(parse[0])
ay2 = float(parse[1])
az2 = float(parse[2])

line1 = file1.readline() #b
parse = line1.split()
bx1 = float(parse[0])
by1 = float(parse[1])
bz1 = float(parse[2])

line2 = file2.readline()
parse = line2.split()
bx2 = float(parse[0])
by2 = float(parse[1])
bz2 = float(parse[2])

line1 = file1.readline() #c
parse = line1.split()
cx1 = float(parse[0])
cy1 = float(parse[1])
cz1 = float(parse[2])

line2 = file2.readline()
parse = line2.split()
cx2 = float(parse[0])
cy2 = float(parse[1])
cz2 = float(parse[2])

line1 = file1.readline() #no.species
line2 = file2.readline()
nspecies = line1.split()

#calc total number of atoms
ntotal = 0
for i in range(len(nspecies)):
	ntotal = ntotal + int(nspecies[i])

line1 = file1.readline() #sel
line2 = file2.readline()
line1 = file1.readline() #cart/direct
line2 = file2.readline()

atomloc1 = []
atomloc2 = []
for i in range (len(nspecies)):
	atomlocsub1 = []
	atomlocsub2 = []
	for j in range (int(nspecies[i])):
		atomlocsub1.append([0]*4)
		atomlocsub2.append([0]*4)
	atomloc1.append(atomlocsub1)
	atomloc2.append(atomlocsub2)

for i in range (len(species)):
	for j in range (int(nspecies[i])):
		line1 = file1.readline()
		if not line1:
			break
		line2 = file2.readline()
		if not line2:
			break

		parse1 = line1.split()
		parse2 = line2.split()
	
		x1 = float(parse1[0])
		y1 = float(parse1[1])
		z1 = float(parse1[2])
	
		x2 = float(parse2[0])
		y2 = float(parse2[1])
		z2 = float(parse2[2])

		atomname = species[i] + str(j+1)

		atomloc1[i][j] = [atomname,x1,y1,z1]
		atomloc2[i][j] = [atomname,x2,y2,z2]

totset = []
for i in range (len(nspecies)):
	j = str(i)
	totset.append(j)

l = 0
parset = []
for cset in xuniqueCombinations(totset,2):
	parset.append(cset)
	l = l+1

for c in range(l):
	a = int(parset[c][0])
	b = int(parset[c][1])

	inidist = []
	for i in range (int(nspecies[a])):
		inidist.append([0]*int(nspecies[b]))

	for i in range (int(nspecies[a])):
		ix = float(atomloc1[a][i][1])
		iy = float(atomloc1[a][i][2])
		iz = float(atomloc1[a][i][3])
		for j in range (int(nspecies[b])):
			jx = float(atomloc1[b][j][1])
			jy = float(atomloc1[b][j][2])
			jz = float(atomloc1[b][j][3])
			vec = [(ix-jx),(iy-jy),(iz-jz)]
			len = length(vec)
			inidist[i][j] = len

	findist = []
	for i in range (int(nspecies[a])):
		findist.append([0]*int(nspecies[b]))

	for i in range (int(nspecies[a])):
		ix = float(atomloc2[a][i][1])
		iy = float(atomloc2[a][i][2])
		iz = float(atomloc2[a][i][3])
		for j in range (int(nspecies[b])):
			jx = float(atomloc2[b][j][1])
			jy = float(atomloc2[b][j][2])
			jz = float(atomloc2[b][j][3])
			vec = [(ix-jx),(iy-jy),(iz-jz)]
			len = length(vec)
			findist[i][j] = len

	print "-------------------------------------------------------------"
	print "  %s - %s  \t ini.distance \t fin.distance \t delta" % (species[a],species[b])
	print "-------------------------------------------------------------"
	for i in range (int(nspecies[a])):
		for j in range (int(nspecies[b])):
			if inidist[i][j] < cutoff:
				delta = float(findist[i][j])-float(inidist[i][j])
				print "%s%i - %s%i \t %5.8f \t %5.8f \t %5.8f" % (species[a],i,species[b],j,inidist[i][j],findist[i][j],delta)

file1.close
file2.close
