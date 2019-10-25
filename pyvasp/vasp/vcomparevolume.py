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

inifile = 'POSCAR'
finfile = 'CONTCAR'

if not os.path.isfile(inifile):
	print(inifile + " is not found.")
	sys.exit()

if not os.path.isfile('CONTCAR'):
	print("\nCONTCAR is not found.")
	finfile = 'POSCAR'

if os.path.isfile('CONTCAR'):
	if not os.path.getsize('CONTCAR'):
		print("\nCONTCAR is empty.")
		finfile = 'POSCAR'

posfile = open(inifile)
line = posfile.readline()
line = posfile.readline()
line = posfile.readline()
parse = line.split()
ax = float(parse[0])
ay = float(parse[1])
az = float(parse[2])
veca = [ax,ay,az]
plena = length(veca)
line = posfile.readline()
parse = line.split()
bx = float(parse[0])
by = float(parse[1])
bz = float(parse[2])
vecb = [bx,by,bz]
plenb = length(vecb)
line = posfile.readline()
parse = line.split()
cx = float(parse[0])
cy = float(parse[1])
cz = float(parse[2])
vecc = [cx,cy,cz]
plenc = length(vecc)
posfile.close()

bcc = crossproduct(vecb, vecc)
pvol = dotproduct(veca, bcc)


confile = open(finfile)
line = confile.readline()
line = confile.readline()
line = confile.readline()
parse = line.split()
ax = float(parse[0])
ay = float(parse[1])
az = float(parse[2])
veca = [ax,ay,az]
clena = length(veca)
line = confile.readline()
parse = line.split()
bx = float(parse[0])
by = float(parse[1])
bz = float(parse[2])
vecb = [bx,by,bz]
clenb = length(vecb)
line = confile.readline()
parse = line.split()
cx = float(parse[0])
cy = float(parse[1])
cz = float(parse[2])
vecc = [cx,cy,cz]
clenc = length(vecc)
confile.close()

bcc = crossproduct(vecb, vecc)
cvol = dotproduct(veca, bcc)

adiff = (clena - plena) / plena * 100
bdiff = (clenb - plenb) / plenb * 100
cdiff = (clenc - plenc) / plenc * 100
voldiff = (cvol - pvol) / pvol * 100

print("              a          b          c        Volume")
print("------------------------------------------------------")
print("POSCAR  : %8.5f   %8.5f   %8.5f   %10.5f" % (plena, plenb, plenc, pvol))
print("CONTCAR : %8.5f   %8.5f   %8.5f   %10.5f" % (clena, clenb, clenc, cvol))
print("DIFF.( ): %8.5f   %8.5f   %8.5f   %10.5f" % (adiff, bdiff, cdiff, voldiff))
