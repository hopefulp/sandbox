#!/opt/applic/epd/bin/python

import sys
import os
import bgf

def is_CO2(myBGF, index):
	"""
is_CO2(myBGF, index):
	check if this atom is contained in water molecule or not and returns the aNo of a water molecule.
	"""
	dest_atom = myBGF.getAtom(index);

	# for C
	if "C" in dest_atom.ffType:
		# it has only two atoms
		if len(dest_atom.CONECT) == 2:
			# ..and they are all Oxygens
			if "O" in myBGF.getAtom(dest_atom.CONECT[0]).ffType and "O" in myBGF.getAtom(dest_atom.CONECT[1]).ffType:
				# ..and their oxygens have no connected atoms
				#if len(myBGF.getAtom(dest_atom.CONECT[0]).CONECT) == 1 and len(myBGF.getAtom(dest_atom.CONECT[1]).CONECT) == 1:
				#	return [dest_atom.aNo, myBGF.getAtom(dest_atom.CONECT[0]).aNo, myBGF.getAtom(dest_atom.CONECT[1]).aNo]
				return [dest_atom.aNo, myBGF.getAtom(dest_atom.CONECT[0]).aNo, myBGF.getAtom(dest_atom.CONECT[1]).aNo]
	# for Oxygen
	elif "O" in dest_atom.ffType:
		alist = dest_atom.CONECT
		if len(alist) == 1:
			if is_CO2(myBGF, myBGF.getAtom(alist[0]).aNo):
				# O H1 H2
				return [myBGF.getAtom(alist[0]).aNo, myBGF.getAtom(alist[0]).CONECT[0], myBGF.getAtom(alist[0]).CONECT[1]]
			return [];
	else:
		return [];


