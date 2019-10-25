#!/opt/applic/epd/bin/python

import sys, re, string, getopt, optparse, math, time, decimal
from os import popen

#-----------------
# move the molecule
# 
#_________________
def movebgf(bgf_file, out_file, x, y, z):
	# read the original bgf
	print(options)

	newcoord = []
	decimal.getcontext().prec = 5

	f_bgf_file = open(bgf_file)
	f_out_file = open(out_file,'w')

	t = time.gmtime()
	f_out_file.write("REMARK moved by moveBGF.py by noische on " + time.asctime(t) + "\n")

	while 1:
		line = f_bgf_file.readline()

		if not line:
			break

		if 'HETATM' or 'ATOM' and not 'FORMAT' in line:
			parse = re.split('\s*',line)
			#parse[6] = str(float(parse[6]) + float(x))
			#wline = '{0:>6} {1:>5} {2:<5} {3:3} {4:<1}{5:>5} {6:>10.5f}{7:>10.5f}{8:>10.5f} {9:<5}{10:3}{11:2} {12:>8.5f}'.format(*item)
			newcoord = [(float(parse[6]) + x), (float(parse[7]) + y), (float(parse[8]) + z)]
			
			wline = '{0:>6}{1:>6} {2:<5} {3:3} {4:<3} {5:3}'.format(parse[0], parse[1], parse[2], parse[3], parse[4], parse[5])
			f_out_file.write(wline)
			fline = '{0:>10.5f}{1:>10.5f}{2:>10.5f} '.format(*newcoord)
			f_out_file.write(fline)
			wline = '{0:<5}{1:3}{2:2} {3:>8.5f}'.format(parse[9], parse[10], parse[11], float(parse[12]))
			f_out_file.write(wline)
			wline = '\n'
			f_out_file.write(wline)
		else:
			f_out_file.write(line)

# if HETATM parse add and copy
# if not HETATM pass

	return 1

# main call

if __name__ == "__main__":

	option = ""; args = ""; bgf_file = ""; mod_file = ""; out_file = ""
	usage = """
	Usage: moveBGF.py -b bgf_file -x coord -y coord -z coord -o out_file
	"""
	
	options, args = getopt.getopt(sys.argv[1:], 'hb:x:y:z:o:', ['help','bgf=','x=','y=','z=','out='])
	for option, value in options:
	        if option in ('-h', '--help'):
	                print usage; sys.exit(0)
	        elif option in ('-b', '--bgf'):
	                bgf_file = value
		elif option in ('-x', '--x'):
			x = float(value)
		elif option in ('-y', '--y'):
			y = float(value)
		elif option in ('-z', '--z'):
			z = float(value)
	        elif option in ('-o', '--out'):
	                out_file = value
	        elif option in (''):
	                print usage; sys.exit(0)
	
	
	movebgf(bgf_file, out_file, x, y, z)
