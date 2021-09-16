#!/usr/bin/env python
import aselite
from sys import argv
import re

if '-h' in argv or len(argv) != 2:
    print('usage: pos2xyz.py POSCAR')
    print()
    exit(1)

fname = argv[1]
atoms = aselite.read_vasp(fname)
if re.match('POSCAR.', fname) or re.match('CONTCAR.', fname):
    fnamelist=fname.split('.')
    outf=fnamelist[1]+'.xyz'
else:
    outf=fname+'.xyz'
print(f"output fname = {outf}")
aselite.write_xyz('%s' % outf, atoms)
