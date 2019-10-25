#!/home/noische/program/python27/bin/python

import sys
import os
import string
import getopt
import nutils as nu
import bgf
import bgftools

version = '110101'

list_aNo = []

def getmolecule(myBGF, atom, list_aNo):
	if atom.aNo not in list_aNo:
		list_aNo.append(atom.aNo)
	else:
		return list_aNo;

	for i in atom.CONECT:
		nextatom = myBGF.getAtom(i)
		getmolecule(myBGF, nextatom, list_aNo)


#-------------------------------------
#
# - Main Function
#
if __name__ == '__main__':

	if len(sys.argv) < 2:
		print(sys.argv[0] + " bgffile aNo")
		sys.exit(1)

	myBGF = bgf.BgfFile(sys.argv[1])
	aNo = int(sys.argv[2])
	atom = myBGF.getAtom(aNo)
	molecule_list = []
	bgftools.getmolecule(myBGF, atom, molecule_list)

	molecule_list.sort()
	print(molecule_list)
