#!/home/jackjack5/epd/bin/python
import math
import os
import sys
import re
import numpy
import scipy.integrate
import getopt
import optparse
import nutils as nu

# by jackjack5 at 20130821

version = "20131101"
R = 0.0;	# radius
sOC = 0.0; 	# sigma_OC parameter
L = 0.0; 	# CNT length
u = 0.0;	# water viscosity
kb = 1.381e-23	# kB
T = 300.0;	# temperature
pi = math.pi;	# pi
forceFile = ""	# force profile (input)
outFilePrefix = ""	# output file prefix ( prefix.out, prefix.dump )
averageStep = 0	# number of steps that will be averaged to calculate AC
sampleStep = 0	# number of steps used to calculate lambda from force AC
averageACStep = 0	# number of steps that will be averaged to calculate lambda from integrated force AC
step = 0	# for how many samples will be extracted


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
	vr = v[averageStep-sampleStep::-1]
	acs = numpy.convolve(v, vr, 'valid')
	acs_av = []	
	n=len(vr)
	for i,a in enumerate(acs):
		temp = a/n
		#temp = a/acs[0]	# for normalization
		acs_av.append(temp)
	return acs_av


def main():
	filename = forceFile
	v = read_value(filename)
	nv = len(v)
	t = numpy.arange(0, nv-1, step)
	result = [];
	result2 = [];
	output = "";

	### calculate constants
	reff = R - 0.5 * sOC
	A = 2 * pi * reff * L
	A *= 1e-20
	C = 1 / kb / T / A

	output += "Version " + str(version) + "\n"
	output += "Parameters\n"
	output += "R" + '\t' + str(R) + '\t' + \
		"r_eff" + '\t' + str(reff) + '\t' + \
		"area" + '\t' + str(A) + '\t' + \
		"coeff" + '\t' + str(C) + '\n' 

	output += str("t") + '\t' + str("int fac") + '\t' + str("l") + '\t' + str("b") + '\n'


	### calculate AC
	print("Reading force file and calculating autocorrelation function..")
	print(str(t))
	for i in t:
		if i + averageStep <= nv:
			temp_v = v[i:i+averageStep]
			temp_acs = ac(temp_v)
			result.append(temp_acs)


	### test for vac
	for index, i in enumerate(result):
		print "** timestep ", t[index]
		a = scipy.integrate.cumtrapz(i, dx=2e-15)	# integration
		a2, _ = nu.meanstdv(a[averageACStep:])
		a2 *= (6.95e-11)**2	# kcal/mol.m -> J/m

		print "\t- integrated force AC=", a2
		l = C * a2
		print "\t- lambda=", l
		b = u / l
		print "\t- slip length (m)=", b
		b *= 1e9
		print "\t- slip length (nm)=", b
		result2.append(b)

		output += str(t[index]) + '\t' + str(a2) + '\t' + str(l) + '\t' + str(b) + '\n'

	
	### write output file
	fname = outFilePrefix + ".a." + str(averageStep) + ".s." + str(sampleStep) + ".n." + str(step) + ".t." + str(averageACStep)
	f_out = open(fname + ".out", 'w')
	f_out.write(str(sys.argv)+"\n")
	f_out.write(output)
	f_out.close()


	### print averaged slip length
	m, _ = nu.meanstdv(result2)
	print "averaged slip length=", m, _;


	### write dump file
	f_dump = open(fname + ".dump", 'w')
	output = "\t";
	# timesteps
	for i in t:
		output += str(i) + "\t"
	output += "\n"
	# result
	for r, row in enumerate(result[0]):
		output += str(r) + "\t"
		for c, col in enumerate(result):
			output += str(result[c][r]) + "\t"
		output += "\n"
	f_dump.write(output)
	f_dump.close()
			
	return 1


if __name__ == "__main__":

	usage = """
	CNT_getVecAC3.py -f force dump file -R radius -r sigma_OC -L CNT length -u viscosity -a averaging steps -s steps for AC -n #timestep for AC -t int range for AC -o prefix 
	"""

	options, args = getopt.getopt(sys.argv[1:], 'hf:R:r:L:u:o:a:s:n:t:', ['help', 'force=', 'R=', 'sigma_OC=', 'length=', 'viscosity=', 'output=', 'average=', 'sample=', 'step=', 'sample2='])
	if len(sys.argv) < 2:
		print(usage)
		sys.exit(0)

	print "Requested options: " + str(options)

	for option, value in options:
	        if option in ('-h', '--help'):
	                print(usage)
			sys.exit(0)
		elif option in ('-f', '--force'):
			forceFile = value
	        elif option in ('-R', '--R'):
			R = float(value)
		elif option in ('-r', '--sigma_OC='):
			sOC = float(value)
		elif option in ('-L', '--length'):
			L = float(value)
		elif option in ('-u', '--viscosity'):
			u = float(value)
		elif option in ('-o', '--prefix'):
			outFilePrefix = value
		elif option in ('-a', '--average'):
			averageStep = int(value)
		elif option in ('-s', '--sample'):
			sampleStep = int(value)
		elif option in ('-n', '--step'):
			step = int(value)
		elif option in ('-t', '--sample2'):
			averageACStep = int(value)

	main()

