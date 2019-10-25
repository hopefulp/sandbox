#!/opt/applic/epd/bin/python

import sys
import os
import numpy

_ = sys.argv[1]
f = open(_)

ll = []
while 1:
	l = f.readline()

	if not l: break;
	if "#" in l: continue;

	l = l.split()

	temp = []
	for i in l:
		temp.append(float(i))
	ll.append(temp)

a = numpy.array(ll)
#print(a)
print numpy.mean(a, axis=0)
