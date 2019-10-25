#!/home/noische/python

import sys
import os

import bgf
import bgftools
import dreiding

usage = "getMass.py bgffile fffile"
if len(sys.argv) < 2:
	print(usage)
	sys.exit(0)

bgf_file = sys.argv[1]
ff_file = sys.argv[2]

myBGF = bgf.BgfFile(bgf_file)
Mw = bgftools.getMass(myBGF, ff_file)

print(bgf_file + " mass: " + str(Mw))
