#!/usr/bin/env python
# -*- coding: utf-8 -*-
###########################################
#
#         COLLECT_DATA 
#
#   script to write reference data into
#   a hdf5 file
#
#    written by JPD (2015)
#
###########################################


import string
import sys
import getopt
import os
import turbomole
import h5py


argv = sys.argv[1:]

def print_help():
    print 'Usage: collect_data -i <reference file>'
    
try:
    opts, args = getopt.getopt(argv, "hpfi:o:")
except getopt.GetoptError:
    print "Unknown option %s requested" % opt
    print_help()
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print_help()
        sys.exit()
    elif opt == '-i':
        reffile = arg 
    else:
        print_help()
        sys.exit(1)

if os.path.isfile('%s' % reffile):
    f = h5py.File('%s' % reffile, 'a')
else:
    print "Reference File not found"
    sys.exit(1)


istruc = string.atoi(os.environ["struc"])
tag = os.environ["tag"]
if os.environ.has_key('cycle'):
    cycle = string.atoi(os.environ['cycle'])
natoms = f['system']['natoms'].value
elements = f['system']['elements'].value
xyz = f['forcematch'][tag]['structures'][istruc,:,:]
tb = turbomole.turbomole(xyz, elements)
tb.cycle = cycle

f['forcematch'][tag]['energies'][istruc] = tb.read_energy()
f['forcematch'][tag]['forces'][istruc,:,:] = tb.read_gradient()

print 'Data successfully written to %s' % reffile
