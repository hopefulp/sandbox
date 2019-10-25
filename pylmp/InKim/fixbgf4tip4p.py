#!/opt/applic/epd/bin/python

import sys
import string

import bgf
import bgftools

myBGF = bgf.BgfFile(sys.argv[1])
fixedBGF = bgf.BgfFile()

for atom in myBGF.a:
	if "OW" in atom.ffType:
		l_watermolecule = atom.CONECT
		fixedBGF.addAtom(atom)	# add OW atom
		for ano in l_watermolecule:
			fixedBGF.addAtom(myBGF.getAtom(ano))	# add two HW atoms
for atom in myBGF.a:
	if atom.ffType != "OW" and atom.ffType != "HW":
		fixedBGF.addAtom(atom)

fixedBGF.saveBGF(sys.argv[1].split(".bgf")[0] + "_mod_tip4p.bgf")

