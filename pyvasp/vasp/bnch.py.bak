#!/usr/bin/python2.7
#heejin
import sys
import os
import re
import subprocess
import operator
import textwrap
import collections

#make force file
if not os.path.isfile('OUTCAR'):
	print "OUTCAR is not found."
	sys.exit()

if len(sys.argv) > 1:
	nstep = int(sys.argv[1])
else:
	nstep = int(100)

command = "grep 'LOOP:  cpu time' OUTCAR | awk '{print $7}' |head -n %d > ~/bnchtemp" % nstep
os.system(command)

homedir = os.path.expanduser('~')
ffile = open(homedir+'/bnchtemp')
print "%5s %10s %10s" % ('step','time','sum')
print "---------------------------"
i = int(0)
sum = float(0)
while 1:
	line = ffile.readline()
	if not line:
		break
	i = i + 1
	sum = sum + float(line)
	print "%5d %10.4f %10.4f" % (i, float(line), sum)

avg = sum / float(i)

print "---------------------------"
print "  Total %10.4f sec" % (sum)
print "Average %10.4f sec/step" % (avg)

ffile.close()
os.system('rm ~/bnchtemp')
