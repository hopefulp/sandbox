#!/usr/bin/python2.7
#heejin
import sys
import os
import re
import math

#usage description
if len(sys.argv)<3:
	print("Usage: [contfile] [chg file]")
	sys.exit()

confilename=sys.argv[1]
chgfilename=sys.argv[2]

# get ZVAL
os.system('grep ZVAL POTCAR | awk \'{print $6}\' > zval.tmp')

# direct to cartesian
poscarfile = open(confilename)
line = poscarfile.readline()
line = poscarfile.readline()
line = poscarfile.readline()
line = poscarfile.readline()
line = poscarfile.readline()
line = poscarfile.readline()
parse = line.split()
if parse[0].isalpha():
	line = poscarfile.readline()
	line = poscarfile.readline()
	line = poscarfile.readline()
	if ('Direct' in line):
		poscarfile.close()
		doos = 'cp ' + confilename + ' CONTCAR.tmp'
		os.system(doos)
		os.system('sed -i 6d CONTCAR.tmp')
		os.system('convasp -cart < CONTCAR.tmp > CONTCAR.cart')
		os.system('rm CONTCAR.tmp')
poscarfile.close()

# get poscar information
poscarfile = open('CONTCAR.cart')

while 1:
	line = poscarfile.readline()
	if not line:
		break
	atomlist = line.split()

	line = poscarfile.readline() #mul
	line = poscarfile.readline() #a
	line = poscarfile.readline() #b
	line = poscarfile.readline() #c

	line = poscarfile.readline()
	nspecies = line.count(' ')
	numlist = line.split()
	natoms = 0
	if numlist[0].isdigit():
		for i in range(nspecies):
			natoms = natoms + int(numlist[i])
	else:
		line = poscarfile.readline()
		nspecies = line.count(' ')
		numlist = line.split()
		for i in range(nspecies):
			natoms = natoms + int(numlist[i])
	break
poscarfile.close()


zval = []
zvalfile = open('zval.tmp')
for i in range(nspecies):
	line = zvalfile.readline()
	zval.append(line)

chgtmpfile = open('CONTCAR.chgtmp', 'w')
chgfile = open(chgfilename)
line = chgfile.readline()
for i in range(nspecies):
	for j in range(int(numlist[i])):
		line = chgfile.readline()
		chglist = line.split()
		chgval = float(zval[i]) - float(chglist[1])
		if chgval > 0:
			chgval2 = '+'+str(chgval)+'\n'
		else:		
			chgval2 = str(chgval)+'\n'
		chgtmpfile.write(chgval2)
chgfile.close()
chgtmpfile.close()

poscarfile = open('CONTCAR.cart')
outfile = open(confilename+'.vo', 'w')
chgfile = open('CONTCAR.chgtmp')

line = poscarfile.readline()
outfile.write(line)
line = poscarfile.readline()
outfile.write(line)
line = poscarfile.readline()
outfile.write(line)
line = poscarfile.readline()
outfile.write(line)
line = poscarfile.readline()
outfile.write(line)
line = poscarfile.readline()
outfile.write(line)
line = poscarfile.readline()
outfile.write(line)
line = poscarfile.readline()
outfile.write(line)

for i in range(natoms):
	line = poscarfile.readline()
	line2 = line.rstrip('\n')
	chgline = chgfile.readline()
	outline = line2 + ' ' + chgline
	outfile.write(outline)

poscarfile.close()
outfile.close()
chgfile.close()

#os.system('convasp -names +2 +5 -2 +1 < CONTCAR.cart > CONTCAR.fc')
os.system('rm zval.tmp CONTCAR.chgtmp CONTCAR.cart')
