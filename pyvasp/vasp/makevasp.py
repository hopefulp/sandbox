#!/usr/bin/python2.7
#heejin
import sys
import os
import re
import math

#enviroments
templetdir = '/qcfs/biduri/templet/vasp'
highk = 50
normalk = 25
highen = 1.5
normalen = 1.3

#usage description
if len(sys.argv)<3:
	print("Usage: [rlx/opt/sp] [high/normal]")
	print(" * rlx-relax atom with fixed cell / opt-cell optimization / sp-single point")
	print(" * high/normal - accuracy level")
	print(" * non-magnetic ordering is default")
	sys.exit()

jobtype=sys.argv[1]
acctype=sys.argv[2]

# direct to cartesian
poscarfile = open('POSCAR')
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
		os.system('cp POSCAR POSCAR.direct')
		os.system('cp POSCAR POSCAR.tmp')
		os.system('sed -i 6d POSCAR.tmp')
		os.system('convasp -cart < POSCAR.tmp > POSCAR')
		os.system('rm POSCAR.tmp')

poscarfile.close()

# get poscar information
poscarfile = open('POSCAR')
def dotproduct(v1, v2):
        return sum((a*b) for a, b in zip(v1, v2))
def length(v):
        return math.sqrt(dotproduct(v, v))

while 1:
	line = poscarfile.readline()
	if not line:
		break
	atomlist = line.split()

	line = poscarfile.readline()

	line = poscarfile.readline()
	x = line.split()
	xvec = [float(x[0]),float(x[1]),float(x[2])]
	xlen = length(xvec)

	line = poscarfile.readline()
	y = line.split()
	yvec = [float(y[0]),float(y[1]),float(y[2])]
	ylen = length(yvec)

	line = poscarfile.readline()
	z = line.split()
	zvec = [float(z[0]),float(z[1]),float(z[2])]
	zlen = length(zvec)

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
#variables : xlen, ylen, zlen, natoms
#for each species: atomlist, numlist

#copy templet files
cptemplet = 'cp '+templetdir+'/INCAR '+templetdir+'/KPOINTS .'
os.system(cptemplet)
if not (os.path.isfile('POTCAR')):
	print("POTCAR file generated\n")
	os.system('potgenerate_pbe POSCAR')
	os.system('mv POSCAR.pot POTCAR')

#modify KPOINTS
if acctype == 'normal':
	kx = int(normalk / int(xlen))
	ky = int(normalk / int(ylen))
	kz = int(normalk / int(zlen))
else:
	kx = int(highk / int(xlen))
	ky = int(highk / int(ylen))
	kz = int(highk / int(zlen))
if kx == 0:
	kx = 1
if ky == 0:
	ky = 1
if kz == 0:
	kz = 1

os.system('sed -i -e \'s|X Y Z|%s %s %s|\' KPOINTS' % (kx, ky, kz))

#get max encut
encuttmp = []
os.system('grep ENMAX POTCAR | awk \'{print $3}\' | cut -c-3 > encuttmp')
tmpfile = open('encuttmp')
while 1:
	line = tmpfile.readline()
	if not line:
		break
	encuttmp.append(line)
encuttmp.sort()
tmpfile.close()
os.system('rm encuttmp')

#modify INCAR
intitle = os.getcwd().split('/')[-1]
inprec = 'normal'
ingga = 'PE'
inismear = '0'
insigma = '0.1'
if acctype == 'high':
	enweight = highen
	inediff = '%sE-5' % str(int(natoms)*1)
	inediffg = '%sE-4' % str(int(natoms)*1)
else:
	enweight = normalen
	inediff = '%sE-5' % str(int(natoms)*2)
	inediffg = '%sE-4' % str(int(natoms)*2)
if int(natoms) < 10:
	inlreal = '.FALSE.'
else:
	inlreal = 'Auto'
if jobtype == 'sp':
	innsw = '0'
	inisif = '2'
elif jobtype == 'opt':
	innsw = '1999'
	inisif = '3'
	inprec = 'high'
else:
	innsw = '1999'
	inisif = '2'
inencut = str(float(encuttmp[-1]) * enweight)
inencut = '520'
inibrion = '2'
inpotim = '0.1'

block = []
magmomlist = ''
for i in range(nspecies):
	block.append(numlist[i]+'*'+'0 ')
	magmomlist = magmomlist+block[i]
os.system('sed -i -e \'s|TITLE|%s|\' INCAR' % intitle)
os.system('sed -i -e \'s|_MAGMOM_|%s|\' INCAR' % magmomlist)
os.system('sed -i -e \'s|_ENCUT_|%s|\' INCAR' % inencut)
os.system('sed -i -e \'s|_PREC_|%s|\' INCAR' % inprec)
os.system('sed -i -e \'s|_ISMEAR_|%s|\' INCAR' % inismear)
os.system('sed -i -e \'s|_SIGMA_|%s|\' INCAR' % insigma)
os.system('sed -i -e \'s|_GGA_|%s|\' INCAR' % ingga)
os.system('sed -i -e \'s|_EDIFF_|%s|\' INCAR' % inediff)
os.system('sed -i -e \'s|_EDIFFG_|%s|\' INCAR' % inediffg)
os.system('sed -i -e \'s|_LREAL_|%s|\' INCAR' % inlreal)
os.system('sed -i -e \'s|_NSW_|%s|\' INCAR' % innsw)
os.system('sed -i -e \'s|_ISIF_|%s|\' INCAR' % inisif)
os.system('sed -i -e \'s|_IBRION_|%s|\' INCAR' % inibrion)
os.system('sed -i -e \'s|_POTIM_|%s|\' INCAR' % inpotim)

print("----- POSCAR information -----")
print("# ATOMS: Total %s atoms" % natoms)
print("LATTICE: %s x %s x %s" % (xlen, ylen, zlen))
print("")
print("----- KPOINT information -----")
print("TYPE   : Monkhorst-Pack")
print("GRID   : %s x %s x %s" % (kx, ky, kz))
print("")
print("----- INCAR  information -----")
print("MAGMOM : %s" % str(magmomlist))
print("ENCUT  : %s" % str(inencut))
print("PREC   : %s" % str(inprec))
print("GGA    : %s" % str(ingga))
print("ISMEAR : %s" % str(inismear))
print("SIGMA  : %s" % str(insigma))
print("EDIFF  : %s" % str(inediff))
print("EDIFFG : %s" % str(inediffg))
print("ISIF   : %s" % str(inisif))
print("IBRION : %s" % str(inibrion))
print("NSW    : %s" % str(innsw))
