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
	return v

def ac(v):
	vr = v[499000::-1]
	acs = numpy.convolve(v, vr, 'valid')
	acs_av = []	
	n=len(vr)
	for i,a in enumerate(acs):
		temp = a/n
		#temp = a/acs[0]	# for normalization
		acs_av.append(temp)
	return acs_av


def main():
	if len(sys.argv)==1:
		print 'Usage : '+sys.argv[0]+' <friction force dump file>'
		return 1
	filename = sys.argv[1]
	v = read_value(filename)
	# assume that force profile is recorded during 2,500,000 steps (5 ns)
	t = [0, 250000, 500000, 750000, 1000000, 1250000, 1500000, 1750000, 2000000]
	result = [];

	### calculate AC
	for i in t:
		temp_v = v[i:i+500000]
		temp_acs = ac(temp_v)
		result.append(temp_acs)

	outfile = open(str(sys.argv[2]) + ".forceAC.1ns.profile", "w")
	output = "\t";

	### timesteps
	for i in t:
		output += str(i) + "\t"
	output += "\n"

	### result
	for r, row in enumerate(result[0]):
		output += str(r) + "\t"
		for c, col in enumerate(result):
			output += str(result[c][r]) + "\t"
		output += "\n"
	outfile.write(output)
	outfile.close()
			
	return 1

main()

