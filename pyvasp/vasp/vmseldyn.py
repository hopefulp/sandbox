#!/usr/bin/python2.7
#heejin
import sys
import os
import re

#usage description
if len(sys.argv)<2:
	posfilename = 'POSCAR'
else:
	posfilename = sys.argv[1]

# write new poscar
posfile = open(posfilename)
newfile = open(posfilename+'.new','w')

stoptk = 0
while 1:
	line = posfile.readline()
	if not line:
		break
	if stoptk is 1:
		break
	newfile.write(line)

	if 'Direct' in line:
		while 1:
			line = posfile.readline()
			if not line:
				break
			if not line.strip():
				stoptk = 1
				break
			parse = line.split()
			outline = "%19.16f %19.16f %19.16f T T T\n" % (float(parse[0]), float(parse[1]), float(parse[2]))
			newfile.write(outline)

posfile.close()
newfile.close()

swfiles = 'mv %s %s.old;mv %s.new %s' % (posfilename, posfilename, posfilename, posfilename)
os.system(swfiles)

chstr = "sed -i -e \'s|Direct|Selective Dynamics\\nDirect|\' %s" % posfilename
os.system(chstr)
