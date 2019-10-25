#!/usr/bin/python2.7
#heejin
import sys
import os
import re

#usage description
if len(sys.argv)<3:
	print("Usage: [file name] [Atom #]")
	sys.exit()

chgfilename=sys.argv[1]
atomnumber=sys.argv[2]

# 
chgfile = open(chgfilename)
tok = 'Atom %s\n' % atomnumber

while 1:
	line = chgfile.readline()
	if not line:
		break

	if tok in line:
		print("%s" % line, end=' ')
		for i in range(60):
			line = chgfile.readline()
			print('%s' % line, end=' ')

chgfile.close()
