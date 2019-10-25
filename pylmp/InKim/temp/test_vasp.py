#!/home/noische/python

import sys
import nutils as nu
from VASP import *

if len(sys.argv) < 2: nu.die(sys.argv[0] + " poscar_file bgf_file")

m = VASP(sys.argv[1])
value = raw_input("Want to create bonds for S? ")
if value: m.find_S_bonds()
m.saveBGF(sys.argv[2])
