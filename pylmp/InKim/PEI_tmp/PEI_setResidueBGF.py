#!/home/noische/program/python27/bin/python

import sys, re, string, getopt, optparse, math, time
from os import popen

option = ""; args = ""; bgf_file = ""; number = ""; out_file = ""
usage = """
Usage: setresidueBGF.py -b bgf_file -n residue_number -o out_file
"""

options, args = getopt.getopt(sys.argv[1:], 'hb:n:o:', ['help','bgf=','number=','out='])
for option, value in options:
        if option in ('-h', '--help'):
                print usage; sys.exit(0)
        elif option in ('-b', '--bgf'):
                bgf_file = value
	elif option in ('-n', '--number'):
		number = value
        elif option in ('-o', '--out'):
                out_file = value
        elif option in (''):
                print usage; sys.exit(0)

#-----------------
# syncronize the residue information between the original BGF file and the modified BGF file by Cerius2
# change RES to PRI, SEC, TER, etc.
#_________________
def setresiduebgf(bgf_file, number, out_file):
	n_atoms = 0

	f_bgf_file = open(bgf_file)
	f_out_file = open(out_file, 'w')

	# count the number of atoms
	while 1:
		line = f_bgf_file.readline()
		if not line:
			break

		if 'HETATM' in line:
			n_atoms += 1
	f_bgf_file.seek(0)

	# write the bgf
	t = time.gmtime()
	f_out_file.write("REMARK Updated by setresidueBGF.py by noische on " + time.asctime(t) + "\n")
	while 1:
		line = f_bgf_file.readline()
		if not line:
			break

		# copy the header
		f_out_file.write(line)

		# very simple parsing method.. just split and get it back again with just updating
		if 'FORMAT ATOM' in line:
			for i in xrange(0, n_atoms):
				line = f_bgf_file.readline()
				parse = re.split('\s*', line)
				parse[5] = str(number) # update the residue

				wline = '{0:>6} {1:>5} {2:<5} {3:3} {4:<1} {5:>5} {6:>9} {7:>9} {8:>9} {9:<6} {10:1} {11:1} {12:>8}'.format(*parse)
				wline += '\n'
				f_out_file.write(wline)
		

	print "setresidueBGF.py: the bgf file " + bgf_file + " is updated with the residue by " + number
	return 1

# main call
setresiduebgf(bgf_file, number, out_file)
