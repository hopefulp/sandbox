#!/opt/applic/epd/bin/python

import sys, re, string, getopt, optparse, math, time
from os import popen

option = ""; args = ""; bgf_file = ""; mod_file = ""; out_file = ""
usage = """
Usage: mergeBGF.py -b bgf1_file -c bgf2_file -o out_file
"""

options, args = getopt.getopt(sys.argv[1:], 'hb:c:o:', ['help','bgf1=','bgf2=','out='])
for option, value in options:
        if option in ('-h', '--help'):
                print usage; sys.exit(0)
        elif option in ('-b', '--bgf1'):
                bgf1_file = value
	elif option in ('-c', '--bgf2'):
		bgf2_file = value
        elif option in ('-o', '--out'):
                out_file = value
        elif option in (''):
                print usage; sys.exit(0)

#-----------------
# merge two bgf file
# 
#_________________
def mergebgf(bgf1_file, bgf2_file, out_file):
	print(options)
# read bgf 1 and bgf 2
	f_bgf1_file = open(bgf1_file)
	f_bgf2_file = open(bgf2_file)
	f_out_file = open(out_file,'w')
	bgf1_atom_data = []; bgf2_atom_data = []; bgf1_conect_data = []; bgf2_conect_data = []
	n_atoms_1 = 0; n_atoms_2 = 0

	while 1:
		line = f_bgf1_file.readline()
		if not line:
			break

		if 'HETATM' in line:
			n_atoms_1 += 1
			parse = re.split('\s*', line)
			bgf1_atom_data.append(parse)
		if 'FORMAT' in line:
			continue
		if 'CONECT' in line:
			parse = re.split('\s*', line)
			parse = parse[:-1]
			bgf1_conect_data.append(parse)

	while 1:
		line = f_bgf2_file.readline()
		if not line:
			break

		if 'HETATM' in line:
			n_atoms_2 += 1
			parse = re.split('\s*', line)
			bgf2_atom_data.append(parse)
		if 'FORMAT' in line:
			continue
		if 'CONECT' in line:
			parse = re.split('\s*', line)
			parse = parse[:-1]
			bgf2_conect_data.append(parse)


# add n_atom_1 to atom id of bgf 2
	#margin = int(math.ceil(n_atoms_1 / 10.0)*10)
	#print(margin)
	margin = n_atoms_1

	for atom in bgf2_atom_data:
		atom[1] = str(int(atom[1]) + margin)
	for conect in bgf2_conect_data:
		n_conect = len(conect) 
		for i in xrange(1, n_conect):
			conect[i] = str(int(conect[i]) + margin)

# merge the file sequentially: 1 -> 2
	f_bgf1_file.seek(0)
	f_bgf2_file.seek(0)

	# header
	while 1:
		line = f_bgf1_file.readline()
		if not line:
			break

		if 'HETATM' in line:
			break

		f_out_file.write(line)

	# atom data of bgf1
	for item in bgf1_atom_data:
		item[6] = float(item[6])
		item[7] = float(item[7])
		item[8] = float(item[8])
		item[12] = float(item[12])
		wline = '{0:>6} {1:>5} {2:<5} {3:3} {4:<1}{5:>5} {6:>10.5f}{7:>10.5f}{8:>10.5f} {9:<5}{10:3}{11:2} {12:>8.5f}'.format(*item)
		wline += '\n'
		f_out_file.write(wline)

	# atom data of bgf2
	for item in bgf2_atom_data:
		item[6] = float(item[6])
		item[7] = float(item[7])
		item[8] = float(item[8])
		item[12] = float(item[12])
		wline = '{0:>6} {1:>5} {2:<5} {3:3} {4:<1}{5:>5} {6:>10.5f}{7:>10.5f}{8:>10.5f} {9:<5}{10:3}{11:2} {12:>8.5f}'.format(*item)
		wline += '\n'
		f_out_file.write(wline)
	
	f_out_file.write('FORMAT CONECT (a6,12i6)\n')

	wline = ""
	for item in bgf1_conect_data:
		for i in xrange(0, len(item)):
			wline += '{0:>6}'.format(item[i])
		wline += '\n'
	f_out_file.write(wline)

	wline = ""
	for item in bgf2_conect_data:
		for i in xrange(0, len(item)):
			wline += '{0:>6}'.format(item[i])
		wline += '\n'
	f_out_file.write(wline)
	f_out_file.write("END\n")
	f_out_file.write("")

	f_out_file.close()	
#return 1

# main call
mergebgf(bgf1_file, bgf2_file, out_file)
