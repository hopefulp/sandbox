##Read some atomic data from util/atomic_data.csv file and store 
# in the appropiate lists and/or dictionaries.
# Data has been obtained from http://edu.kde.org/kalzium (version 2.4.00)
import os
import csv
from collections import OrderedDict
import string
#import re

#conv. factor
pm2bohr = 52.917721092
pm2ang  = 0.01

#Read in data and convert to a.u. units
numbers_in = []
labels_in = []
mass_in = [] # a.u.
cov_radii_in = [] # picometers
vdw_radii_in = [] # picometers
with open(os.path.join(os.environ['PDLPDIR'],'py/atomic_data.csv'), 'r') as symb:
	reader = csv.reader(symb, delimiter=',', skipinitialspace=True)
	reader.next()
	for row in reader:
		numbers_in.append(int(row[1]))
		labels_in.append(row[2].lower())
		mass_in.append(row[3].split()[0].replace(",", "."))
		cov_radii_in.append(float(row[5].split()[0]) * pm2ang)
		vdw_radii_in.append(float(row[6].split()[0]) * pm2ang)

def print_dict(dictionary):
	for key, value in dictionary.iteritems():
	    print "%s: %s" % (key, value)

#Labels 
labels_in.insert(0,'xx')
labels = labels_in
# for compatibility with weaver
elems = labels

#Numbers
numbers_in.insert(0,0)
numbers = OrderedDict(zip(labels, numbers_in))

#Masses
mass_in.insert(0,0.0)
mass = OrderedDict(zip(labels, mass_in))

#Covalent radii
cov_radii_in.insert(0,0.3)
cov_radii = OrderedDict(zip(labels, cov_radii_in))

#Covalent radii
vdw_radii_in.insert(0,0.3)
vdw_radii = OrderedDict(zip(labels, vdw_radii_in))

#Just for debugging
#print_dict(numbers)

def call(elem, info):
    infos = {'cov_radii':cov_radii, 'mass':mass, 'number':numbers, 'vdw_radii':vdw_radii}
    elem = string.lower(string.split(elem)[0])
    return infos[info][elem]
