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

if len(sys.argv) < 2:
	print("Usage: [file1] [file2]")
	sys.exit()

filename1 = sys.argv[1]
filename2 = sys.argv[2]

if not os.path.isfile(filename1):
	print(inifile + " is not found.")
	sys.exit()

if not os.path.isfile(filename2):
	print(finfile + " is not found.")
	sys.exit()

file1 = open(filename1)
file2 = open(filename2)

tmpname1 = filename1 + ".tmp"
tmpname2 = filename2 + ".tmp"
tmpfile1 = open(tmpname1,'w')
tmpfile2 = open(tmpname2,'w')

line1 = file1.readline() #species
line2 = file2.readline()
tmpfile1.write(line1)
tmpfile2.write(line2)
species = line1.split()

line1 = file1.readline() #scale
line2 = file2.readline()
tmpfile1.write(line1)
tmpfile2.write(line2)

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

aax = (ax1+ax2)/2 #average
aay = (ay1+ay2)/2
aaz = (az1+az2)/2
abx = (bx1+bx2)/2
aby = (by1+by2)/2
abz = (bz1+bz2)/2
acx = (cx1+cx2)/2
acy = (cy1+cy2)/2
acz = (cz1+cz2)/2

veca = str(aax) + " " + str(aay) + " " + str(aaz) + "\n"
vecb = str(abx) + " " + str(aby) + " " + str(abz) + "\n"
vecc = str(acx) + " " + str(acy) + " " + str(acz) + "\n"

tmpfile1.write(veca)
tmpfile2.write(veca)
tmpfile1.write(vecb)
tmpfile2.write(vecb)
tmpfile1.write(vecc)
tmpfile2.write(vecc)

line1 = file1.readline() #name.species
line2 = file2.readline()

line1 = file1.readline() #no.species
line2 = file2.readline()
tmpfile1.write(line1)
tmpfile2.write(line2)
nspecies = line1.split()

line1 = file1.readline() #sel
line2 = file2.readline()
tmpfile1.write(line1)
tmpfile2.write(line2)
line1 = file1.readline() #cart/direct
line2 = file2.readline()
tmpfile1.write(line1)
tmpfile2.write(line2)

while 1:
	line1 = file1.readline()
	line2 = file2.readline()
	if not line1:
		break
	if not line2:
		break
	tmpfile1.write(line1)
	tmpfile2.write(line2)

file1.close
file2.close
tmpfile1.close
tmpfile2.close

#convert
cname1 = tmpname1 + '.cart'
cname2 = tmpname2 + '.cart'
com1 = 'convasp -cart < ' + tmpname1 + ' > ' + cname1
com2 = 'convasp -cart < ' + tmpname2 + ' > ' + cname2
os.system(com1)
os.system(com2)

# new...
file1 = open(cname1)
file2 = open(cname2)

line1 = file1.readline() #species
line2 = file2.readline()
line1 = file1.readline() #scale
line2 = file2.readline()
line1 = file1.readline() #a
line2 = file2.readline()
line1 = file1.readline() #b
line2 = file2.readline()
line1 = file1.readline() #c
line2 = file2.readline()
line1 = file1.readline() #no.species
line2 = file2.readline()
line1 = file1.readline() #sel
line2 = file2.readline()
line1 = file1.readline() #cart/direct
line2 = file2.readline()

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

		if (x2-x1)>5:
			x2 = x2-aax
			y2 = y2-aay
			z2 = z2-aaz
		if (x2-x1)<-5:
			x2 = x2+aax
			y2 = y2+aay
			z2 = z2+aaz
		if (y2-y1)>5:
			x2 = x2-abx
			y2 = y2-aby
			z2 = z2-abz
		if (y2-y1)<-5:
			x2 = x2+abx
			y2 = y2+aby
			z2 = z2+abz
		if (z2-z1)>5:
			x2 = x2-acx
			y2 = y2-acy
			z2 = z2-acz
		if (z2-z1)<-5:
			x2 = x2+acx
			y2 = y2+acy
			z2 = z2+acz

		vec = [(x2-x1),(y2-y1),(z2-z1)]
		len = length(vec)
		atname = species[i]+str(j+1)
	#	print "%-5s %8.5f from( %8.5f , %8.5f , %8.5f ) to( %8.5f , %8.5f , %8.5f )" % (atname,len,x1,y1,z1,x2,y2,z2)
		print("%-5s %8.5f" % (atname,len))

file1.close
file2.close

comm = 'rm ' + tmpname1 + ' ' + tmpname2 + ' ' + cname1 + ' ' + cname2
#os.system(comm)
