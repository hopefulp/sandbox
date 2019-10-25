#!/usr/bin/python2.7
#heejin
import sys
import os
import re

#usage description
if len(sys.argv)<3:
	print("Usage: [POSCAR file] [del list file]")
	sys.exit()

posfilename=sys.argv[1]
delfilename=sys.argv[2]

# get del list
delfile = open(delfilename)
dellist = []
while 1:
	line = delfile.readline()
	if not line:
		break
	dellist.append(int(line))

# write new poscar
posfile = open(posfilename)
newfile = open(posfilename+'.new','w')
while 1:
	line = posfile.readline()
	if not line:
		break
	newfile.write(line)

	if 'Direct' in line:
		i = 0
		while 1:
			line = posfile.readline()
			i = i + 1
			if not line:
				break
			if i in dellist:
				print("Line %d removed" % i)
			else:
				newfile.write(line)
swfiles = 'mv %s %s.old;mv %s.new %s' % (posfilename, posfilename, posfilename, posfilename)
os.system(swfiles)
