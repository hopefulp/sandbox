#!/usr/bin/python2.7
#heejin
import sys
import os
import re
import math
import subprocess
import operator

#usage description
if len(sys.argv) < 4:
	print "Usage: [POSCAR file] [Atom name] [# remaining] [sub. name]"
	sys.exit()

posfilename=sys.argv[1]
targetname=sys.argv[2]
nsubstr=sys.argv[3]
subname=sys.argv[4]

if not os.path.isfile(posfilename):
	print posfilename + " is not found."
	sys.exit()
if not os.path.isfile('INCAR'):
	print "INCAR is not found."
	sys.exit()
if not os.path.isfile('KPOINTS'):
	print "KPOINTS is not found."
	sys.exit()
if nsubstr.isdigit():
	nsubatom = nsubstr
else:
	nsubatom = 0
	print "# of remaining atoms is wrong."
	sys.exit()

#read atom species from POTCAR
if not os.path.isfile('POTCAR'):
	print "POTCAR is not found."
	sys.exit()
else:
	species = []
	for line in open('POTCAR'):
		if "TITEL" in line:
			sptemp = line.split()[3]
			species.append(re.split('_',sptemp)[0])

print species

#read poscar
poscarfile = open(posfilename)
for i in range (5):
	poscarfile.readline()
line = poscarfile.readline()
if line.split()[0].isalpha():
	line = poscarfile.readline()
	nspecies = line.split()
else:
	nspecies = line.split()
poscarfile.close()

if len(nspecies) != len(species):
	print posfilename +" and POTCAR are inconsistent."
	sys.exit()

#calc number of atoms
ntotal = 0
ntarget = 0
for i in range(len(nspecies)):
	if species[i] in targetname:
		ntarget = int(nspecies[i])
	ntotal = ntotal + int(nspecies[i])
if ntarget is 0:
	print targetname +" is not in "+ posfilename
	sys.exit()

#make combination list
def xuniqueCombinations(items, n):
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for cc in xuniqueCombinations(items[i+1:],n-1):
                yield [items[i]]+cc\

combfile = open('combination.tmp','w')
totset = []
for j in range (ntarget):
	k = str(j + 1)
	totset.append(k)

for c in xuniqueCombinations(totset,int(nsubatom)):
	combfile.write(' '.join(c))
	combfile.write('\n')
combfile.close()

##-build poscar
poscarfile = open(posfilename)
outfile = open(posfilename+'.head', 'w')
outfile2 = open(posfilename+'.target', 'w')

line = poscarfile.readline()
line2 = line.rstrip('\n')
outline = line2 + ' ' + subname
outfile.write(outline+'\n')

for i in range(4):
	line = poscarfile.readline()
	outfile.write(line)

line = poscarfile.readline()
line2 = line.rstrip('\n')
outline = line2 + ' ' + subname
outfile.write(outline+'\n')

line = poscarfile.readline()
outline = ''
nsubatom2 = int(ntarget) - int(nsubatom)
for i in range (len(nspecies)):
	if species[i] in targetname:
		outline = outline + ' ' + nsubatom + ' ' + str(nsubatom2)
	else:
		outline = outline + ' ' + str(nspecies[i])
outfile.write(outline+'\n')

line = poscarfile.readline()
if len(line.split()) is 2: #selective dynamics
	outfile.write(line)
	line = poscarfile.readline()
	outfile.write(line)
else:
	outfile.write(line)

for i in range((ntotal-ntarget)):
	line = poscarfile.readline()
	line2 = line.rstrip('\n')
	outline = line2
	outfile.write(outline+'\n')

for i in range(ntarget):
	line = poscarfile.readline()
	line2 = line.rstrip('\n')
	outline = line2
	outfile2.write(outline+'\n')

poscarfile.close()
outfile.close()
outfile2.close()

#make configuration
combfile = open('combination.tmp')
targetfile = open(posfilename + '.target')

i = 0
while 1:
	line = combfile.readline()
	i = i + 1
	if not line:
		break
	occlist = line.split()
	print "%d - %s\n" % (i, occlist)
	
	targetfile.close()
	targetfile = open(posfilename + '.target')
	newfilename = posfilename + '_' + str(i)
	newfilenameb = posfilename + '_' + str(i) + 'b'
	newfile = open(newfilename,'w')
	newfileb = open(newfilenameb,'w')
	headfile = open(posfilename + '.head')
	while 1:
		line = headfile.readline()
		if not line:
			break
		newfile.write(line)
	headfile.close()

	k = 0
	while 1:
		line = targetfile.readline()
		if not line:
			break
		k = k + 1
		if str(k) in occlist:
			newfile.write(line)
		else:
			newfileb.write(line)
	newfile.close()
	newfileb.close()

	#append 2nd atoms
	app = 'cat ' + newfilenameb + ' >> ' + newfilename
	os.system(app)

targetfile.close()
combfile.close()


#make directories
tmpfile = open('combination.tmp')
i = 0
while 1:
	line = tmpfile.readline()
	if not line:
		break
	parse = '_'.join(line.split())
	dname = 'con' + str(i) + '_' + parse
	mkdr = "mkdir %s" % dname
	os.system(mkdr)

	pname = posfilename + '_' + str(i+1)
	pname2 = 'INCAR KPOINTS POTCAR'
	tname = dname + '/POSCAR'

	wr = "cp %s %s" % (pname, tname)
	wr2 = "cp %s %s" % (pname2, dname)

	os.system(wr)
	os.system(wr2)
	i = i + 1
tmpfile.close()

#delete temporary files
dd1 = 'rm ' + posfilename + '_*'
os.system(dd1)
os.system('rm *.head *.target')
