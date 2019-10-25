#!/usr/bin/python2.7
#heejin
import os
import sys
import re

#usage description
if len(sys.argv)<3:
	print("Usage: [prefix.in] [prefix.out]")
	print("ex: afmn3-dos-DOS fehse")
	exit()

prefixi=sys.argv[1]
prefixo=sys.argv[2]

# get poscar information
#for each species: atomlist, numlist
#total number of atoms : natoms
poscarfile = open('CONTCAR')
while 1:
	line = poscarfile.readline()
	if not line:
		break
	atomlist = line.split()

	line = poscarfile.readline()
	line = poscarfile.readline()
	line = poscarfile.readline()
	line = poscarfile.readline()

	line = poscarfile.readline()
	numlist = line.split()
	nspecies = len(numlist)
	natoms = 0

	if numlist[0].isdigit():
		for i in range(nspecies):
			natoms = natoms + int(numlist[i])
	else:
		line = poscarfile.readline()
		numlist = line.split()
		nspecies = len(numlist)
		for i in range(nspecies):
			natoms = natoms + int(numlist[i])
	break
poscarfile.close()

#
inin = 1
nelem = []
dolist = []
for i in range(nspecies):
	elem = ''
	for k in range(inin, int(numlist[i])+inin):
		elem = str(elem) + ' ' + str(k)
	inin = inin + int(numlist[i])
	nelem.append(elem)

	doi = 'sumdos2.py ' + prefixi + str(nelem[i]) + ' ' + prefixo + '-' + str(atomlist[i])
	print('Making PDOS for ' + str(atomlist[i]))
	os.system(doi)

for k in range(1, int(natoms)+1):
	elem = str(elem) + ' ' + str(k)
doi = 'sumdos2.py ' + prefixi + elem + ' ' + prefixo + '-Total'
print('Making Total DOS')
os.system(doi)
