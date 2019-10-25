#!/home/jackjack5/epd/bin/python
import math
import os
import sys
import re
import random
import numpy
import time
# by jackjack5 at 20130821

def read_value(filename):
	f=open(filename)
	fsp = f.readlines()
	v = []
	for l in fsp:
		if '#' in l:
			continue;

		lsp = l.split()

		v.append(float(lsp[-1]))
	#print(len(v))
	return v

def ac(v):
	vr = v[499000::-1]
	acs = numpy.convolve(v, vr, 'valid')
	acs_av = []	
	n=len(vr)
	for i,a in enumerate(acs):
		temp = a/n
		#temp = a/acs[0]
		acs_av.append(temp)
		print(str(i) + "\t" + str(temp))
	return acs_av


def main():
	if len(sys.argv)==1:
		print 'Usage : '+sys.argv[0]+' <friction force dump file>'
		return 1
	filename = sys.argv[1]
	v = read_value(filename)
	acs = ac(v)	
	return 1

main()

