#!/home/noische/python

import sys
import os
import getopt
from vasp2bgf import *


def makeBonds():
	pass;

def main():
	pass;

if __name__ == "__main__":
	poscar_file = ""; out_file = ""; 

	options, args = getopt.getopt(sys.argv[1:], 'hp:o:', ['help', 'poscar=', 'out='])

	if len(sys.argv) < 2:
		print(usage)
		sys.exit(0)

	print "Requested options: " + str(options)

	for option, value in options:
	        if option in ('-h', '--help'):
	                print(usage)
			sys.exit(0)
	        elif option in ('-p', '--poscar'):
	                poscar_file = value
	        elif option in ('-o', '--out'):
	                out_file = value
	        elif option == NULL:
			print(usage)
			sys.exit(0)
	
	# main call
	main(poscar_file, out_file)
