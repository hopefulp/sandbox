#!/usr/bin/python2.7
#heejin
import sys
import os
import re

currdir = os.path.basename(os.getcwd())
dochg = '/qcfs/biduri/scripts/vasp/convasp -chgint CHGCAR > '+currdir+'.chg &'

chgcarfile = open('CHGCAR')
line = chgcarfile.readline()
line = chgcarfile.readline()
line = chgcarfile.readline()
line = chgcarfile.readline()
line = chgcarfile.readline()
line = chgcarfile.readline()
parse = line.split()
chgcarfile.close()
if parse[0].isalpha():
	os.system('sed -i 6d CHGCAR')
os.system(dochg)
