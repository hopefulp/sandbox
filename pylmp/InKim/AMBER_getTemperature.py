#!/opt/applic/epd/bin/python

import os
import sys
import string
import re

if len(sys.argv) < 3:
	print("usage: AMBER_getDensity.py mden_file out_file")
	sys.exit(0)

if os.path.exists(sys.argv[2]):
	nu.die("File for output already exists.")

mden = sys.argv[1]
f_mden = open(mden)
out = sys.argv[2]
f_out = open(out, 'w')
step = 0; temperature = 0;

while 1:
	line = f_mden.readline()
	#parse = line.split(' ')
	parse = re.split('\s*', line)

	if not line:
		break

	if 'L0' in parse:
		if "step" in parse[1]:
			step = parse[1]
		else:
			step = int(parse[1])
	if 'L1' in parse:
		if "Temp" in parse[1]:
			temperature = parse[1]
		else:
			temperature = float(parse[1])

	if step != 0 and temperature != 0:
		output = str(step) + "\t" + str(temperature) + "\n"
		#print(output)
		f_out.write(output)
		step = 0; temperature = 0;

### end of function
