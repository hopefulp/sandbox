#!/home/noische/python

import sys
from CNT import *

m=Nanotube(sys.argv[1])
m.make_pbc()
m.make_infinite()
