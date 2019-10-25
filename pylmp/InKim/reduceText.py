#!/opt/applic/epd/bin/python

import os, sys
import nutils as nu

"""
This script read every first column and reduce the line size 
"""
def main():
	usage = """This script averages the data which has rows in 'step number', usually from LAMMPS profile.
	Usage: reduceText.py input output Nevery Nrepeat Nfreq
	With (Nevery, Nrepeat, Nfreq) == (2, 6, 100), the following will be averaged:
		[90, 92, 94, 96, 98, 100] --> 100
		[190, 192, 194, 196, 198, 200] --> 200
		...
	"""
	if len(sys.argv) < 2:
		print("Usage: " + str(sys.argv[0]) + " input output Nevery Nrepeat Nfreq")
		sys.exit(0)

	infile = sys.argv[1];
	outfile = sys.argv[2];
	n_every = int(sys.argv[3]);
	n_repeat = int(sys.argv[4]);
	n_freq = int(sys.argv[5]);

	ifs = open(infile);
	ofs = open(outfile, 'w');

	n_count = 0;	# main counter
	output = "";
	temp = []; 

	while 1:

		line = ifs.readline()
		if not line:
			break;

		if line[0] == "#":
			continue;

		parse = line.split();
		step = int(parse[0])
		value = float(parse[1])

		n_count += 1;

		if n_count <= n_freq:
			#temp.append(step);	#test
			temp.append(value);
		else:
			#temp.append(step);	#test
			#print(temp)		#test
			#print(temp[::n_every][-n_repeat:])	#test
			temp.append(value);
			m, _ = nu.meanstdv(temp[::n_every][-n_repeat:])		# average
			if _ != 0.0:
				output += str(parse[0]) + '\t' + str(m) + '\t' + str(_) + '\n'
			else:
				output += str(parse[0]) + '\t' + str(m) + '\n'


			# initialize
			temp = [0];
			n_count = 1;

	ofs.write(output)
	print("File is written in " + outfile + " ...")

if __name__ == "__main__":
	main()

