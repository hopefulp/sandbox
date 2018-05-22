#!/usr/bin/python2.7
#heejin
import sys
import os
import re
import math
import subprocess
import operator

#usage description
if len(sys.argv) < 5:
	print "Usage: [POSCAR file] [chg file] [Atom name] [# Remaining] [# Max Conf.]"
	sys.exit()

posfilename=sys.argv[1]
chgfilename=sys.argv[2]
targetname=sys.argv[3]
nrmngstr=sys.argv[4]
enmaxstr=sys.argv[5]

if not os.path.isfile(posfilename):
	print posfilename + " is not found."
	sys.exit()
if not os.path.isfile(chgfilename):
	print chgfilename + " is not found."
	sys.exit()
if not os.path.isfile('INCAR'):
	print "INCAR is not found."
	sys.exit()
if not os.path.isfile('KPOINTS'):
	print "KPOINTS is not found."
	sys.exit()
if nrmngstr.isdigit():
	nremaining = nrmngstr
else:
	nremaining = 0
	print "# of remaining atoms is wrong."
	sys.exit()
if enmaxstr.isdigit():
	enmax = int(enmaxstr)
else:
	enmax = 20
	print "# of configurations is wrong."
	sys.exit()

#read atom species from POTCAR
if not os.path.isfile('POTCAR'):
	print "POTCAR is not found."
	sys.exit()
else:
	species = []
	zspecies = []
	for line in open('POTCAR'):
		if "TITEL" in line:
			sptemp = line.split()[3]
			species.append(re.split('_',sptemp)[0])
		if "ZVAL" in line:
			zspecies.append(line.split()[5])

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

for c in xuniqueCombinations(totset,int(nremaining)):
	combfile.write(' '.join(c))
	combfile.write('\n')
combfile.close()

#make poscar with charge
def del6line(file):
	poscarfile.close()
	del6line = 'sed -i 6d ' + file
	os.system(del6line)

def con2cart(file):
	cppos = 'cp ' + file + ' postemp'
	convt = 'convasp -cart < postemp >' + file
	os.system(cppos)
	os.system(convt)
	os.system('rm postemp')

cpposfile = 'cp '+ posfilename +' '+ posfilename+'.org'
os.system(cpposfile)

poscarfile = open(posfilename)
for i in range (5):
	poscarfile.readline()
line = poscarfile.readline()

  ##-convert to cartesian
if line.split()[0].isalpha(): #vesta format
	line = poscarfile.readline()
	line = poscarfile.readline()
	if len(line.split()) is 2: #selective dynamics
		line = poscarfile.readline()
		if 'D' in line: #direct
			del6line(posfilename)
			con2cart(posfilename)
		else: #cart
			del6line(posfilename)
	else: #not selective dynamics
		if 'D' in line: #direct
			del6line(posfilename)
			con2cart(posfilename)
		else: #cart
			del6line(posfilename)
else: #not vesta format
	line = poscarfile.readline()
	if len(line.split()) is 2: #selective dynamics
		line = poscarfile.readline()
		if 'D' in line: #direct
			con2cart(posfilename)
	else: #not selective dynamics
		if 'D' in line: #direct
			con2cart(posfilename)

  ##-build temporary charge info file
chgtmpfile = open('charge.tmp', 'w')
chgfile = open(chgfilename)
line = chgfile.readline()
for i in range(len(nspecies)):
	for j in range(int(nspecies[i])):
		line = chgfile.readline()
		chglist = line.split()
		chgval = float(zspecies[i]) - float(chglist[1])
		if chgval > 0:
			chgval2 = '+'+str(chgval)+'\n'
		else:		
			chgval2 = str(chgval)+'\n'
		chgtmpfile.write(chgval2)
chgfile.close()
chgtmpfile.close()

  ##-build poscar
poscarfile = open(posfilename)
outfile = open(posfilename+'.head', 'w')
outfile2 = open(posfilename+'.target', 'w')
chgfile = open('charge.tmp')

for i in range(5):
	line = poscarfile.readline()
	outfile.write(line)
line = poscarfile.readline()

outline = ''
for i in range (len(nspecies)):
	if species[i] in targetname:
		outline = outline + ' ' + nremaining
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
	chgline = chgfile.readline()
	outline = line2 + ' ' + chgline
	outfile.write(outline)

for i in range(ntarget):
	line = poscarfile.readline()
	line2 = line.rstrip('\n')
	chgline = chgfile.readline()
	outline = line2 + ' ' + chgline
	outfile2.write(outline)

poscarfile.close()
outfile.close()
outfile2.close()
chgfile.close()

#make configuration
combfile = open('combination.tmp')
targetfile = open(posfilename + '.target')

i = 0
outfile = open('ewaldsum.out','w')
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
	newfile = open(newfilename,'w')
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
			pass
	newfile.close()

	ew = "convasp -ewald < %s |tail -n 1" % newfilename
	pp = subprocess.Popen(ew, shell=True, stdout=subprocess.PIPE)
	stdout = pp.stdout
	outline = str(i) + ' ' + stdout.readline()
	outfile.write(outline)
targetfile.close()
combfile.close()
outfile.close()

#pick the energy values
tmpfile = open('ewaldsum.out')
tmptarget = open('ewaldsum.out2','w')
while 1:
	line = tmpfile.readline()
	if not line:
		break
	parse = line.split()
	outline = "%10.6f %s\n" % (float(parse[4]), parse[0])
	tmptarget.write(outline)
tmpfile.close()
tmptarget.close()

#sort ewald-sum
tmpfile = open('ewaldsum.out2')
tmptarget = open('ewaldsum.out3','w')
lines = tmpfile.readlines()
lines.sort(reverse=True)
tmptarget.writelines(lines)
tmpfile.close()
tmptarget.close()

#pick enmax lowest-E configuration
tmpfile = open('ewaldsum.out3')
tmptarget = open('ewaldsum.out4','w')
tempenergy = 'abcdegf'
i = 0
while i < enmax:
	line = tmpfile.readline()
	if not line:
		break
	parse = line.split()

	i = i + 1
	tempenergy = parse[0]
	outline = "%s\n" % parse[1]
	tmptarget.write(outline)
tmpfile.close()
tmptarget.close()

#make directories
tmpfile = open('ewaldsum.out4')
i = 0
while 1:
	line = tmpfile.readline()
	if not line:
		break
	parse = line.split()
	dname = 'con' + str(i) + '_' + parse[0]
	mkdr = "mkdir %s" % dname
	os.system(mkdr)

	pname = posfilename + '_' + parse[0]
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
os.system('rm ewaldsum.out2 ewaldsum.out3 ewaldsum.out4 charge.tmp combination.tmp *.head *.target')
