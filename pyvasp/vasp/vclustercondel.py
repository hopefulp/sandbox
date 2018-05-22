#!/usr/bin/python2.7
#heejin
import sys
import os
import re

contfile = open('CONTCAR')
line = contfile.readline()
line = contfile.readline()
line = contfile.readline()
line = contfile.readline()
line = contfile.readline()
line = contfile.readline()
parse = line.split()
contfile.close()
if parse[0].isalpha():
	os.system('sed -i 6d CONTCAR')
