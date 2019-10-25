#!/home/noische/program/python27/bin/python
"""
example.py
Contains some examples dealing with bgffile and bgfatom classes.
"""

import sys
if '/home/noische/script' not in sys.path:
	sys.path.append('/home/noische/script')
from bgf import *
import bgftools
from updateBGF import *
from centerBGF import *

#myBGF = BgfFile("600_E001_solv.bgf")
#print(bgftools.is_water(myBGF, 115))
#print(bgftools.is_water(myBGF, 126))
#print(bgftools.is_water(myBGF, 108))
#print(len(bgftools.listOxygenAtoms(myBGF)))
#print(len(bgftools.listWaterAtoms(myBGF)))
#bgftools.deleteWaterAtoms(myBGF, 115)
#myBGF.saveBGF("delwatertest.bgf")
#print(len(bgftools.listSoluteAtoms(myBGF)))
#bgftools.removeBadContacts(myBGF, 5.0)

#myBGF.saveBGF("delwatertest.bgf")

#trj_file = str(sys.argv[2])
#out_file = str(sys.argv[3])
#updatebgf(bgf_file, trj_file, out_file)
bgf_file = str(sys.argv[1])
out_file = str(sys.argv[2])
myBGF = BgfFile(bgf_file)
myBGF2 = bgftools.replicateCell(myBGF, (3, 3, 3), True)
myBGF.saveBGF(out_file)
