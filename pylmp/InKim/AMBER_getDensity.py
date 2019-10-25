#!/opt/applic/epd/bin/python

import os
import sys
import string
import re
import time
import nutils as nu

version = 'kdft121011'

if len(sys.argv) < 3:
	print("usage: AMBER_getDensity.py mden_file top_file out_file")
	sys.exit(0)

if os.path.exists(sys.argv[3]):
	nu.die("File for output already exists.")

mden = sys.argv[1]
top = sys.argv[2]
out = sys.argv[3]

f_mden = open(mden)
f_top = open(top)
f_out = open(out, 'w')

f_out.write("Requested options: " + str(sys.argv) + "\n")
f_out.write("Job started at " + time.asctime(time.gmtime()) + " on " + os.environ["HOSTNAME"] + " by " + os.environ["USER"] + " at " + os.environ["PWD"] + "\n")
step = 0; density = 0; volume = 0; temp = 0;

### parse top
parse_mass = [];
while 1:
	line = f_top.readline()

	if not line:
		break;

	if 'FLAG MASS' in line:
		while 1:
			line2 = f_top.readline()
			parse = re.split("\s*", line2)

			if 'FLAG' in line2:
				break;
			else:
				parse_mass.append(parse)

l_mass = nu.flatten(parse_mass)

# add mass
mass = 0;
for i in l_mass:
	try:
		mass += float(i)
	except:
		pass;

#print(mass)


### parse mden
while 1:
	line = f_mden.readline()
	parse = re.split('\s*', line)

	if not line:
		break;

	if 'L0' in parse:
		if "step" in parse[1]:
			step = parse[1]
		else:
			step = int(parse[1])

	if 'L1' in parse:
		if "Temp" in parse[1]:
			temp = parse[1]
		else:
			temp = float(parse[1])

	if 'L3' in parse:
		if "vol" in parse[1]:
			volume = parse[1]
		else:
			volume = float(parse[1])

	if 'L9' in parse:
		if "Dens" in parse[4]:
			density = parse[4]
		else:
			density = float(parse[4])

	if temp != 0 and step != 0 and density != 0:
		output = str(step) + "\t" + str(temp) + "\t" + str(volume) + "\t" + str(mass) + "\t" + str(density) + "\n"
		#print(output)
		f_out.write(output)
		step = 0; density = 0; volume = 0; temp = 0;

### end of function
