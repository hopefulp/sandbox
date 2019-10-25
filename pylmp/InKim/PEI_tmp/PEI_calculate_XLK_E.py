#!/usr/bin/python
"""
run singlepoint calculation with these:

compute pbl all property/local btype batom1 batom2
compute pal all property/local atype aatom1 aatom2 aatom3
compute pdl all property/local dtype datom1 datom2 datom3 datom4
compute bl all bond/local dist eng
compute al all angle/local theta eng
dump    1 all local 1 bond.profile index c_pbl[1] c_pbl[2] c_pbl[3] c_bl[1] c_bl[2]
dump    2 all local 1 angle.profile index c_pal[1] c_pal[2] c_pal[3] c_pal[4] c_al[1] c_al[2]
run 0

"""

import sys
import string
import bgf

def equal_split(list, n):
	L = []
	while list != []:
		L.append(list[:n])
		list = list[n:]

	return L

b = bgf.BgfFile(sys.argv[1])
l = [];
c1 = ""; c5 = ""; n = ""; m = "";
xlnktype = dict()

for atom in b.a:
	if atom.rName == "XLK":
		l.append(atom.aNo)

	if "XLK" in atom.rName and "C1" in atom.aName:
		c1 = atom
		c5 = b.getAtom(b.a2i[c1.aNo + 3])

		for number in c1.CONECT:
			temp = b.getAtom(number)
			if "N" in temp.aName:
				if bgf.is_bonded(temp, c1):
					n = temp

		for number in c5.CONECT:
			temp = b.getAtom(number)
			if "N" in temp.aName:
				if bgf.is_bonded(temp, c5):
					m = temp

		if n.rNo != m.rNo:
			dist = bgf.distance(n, m)
			xlnktype[c1.aNo] = "inter"
		else:
			dist = bgf.distance(n, m)
			xlnktype[c1.aNo] = "intra"

L = equal_split(l, 10)

# energy
E = [0 for i in L];

bondinfo = open(sys.argv[2])
angleinfo = open(sys.argv[3])

#bondinfo
while 1:
	line = bondinfo.readline()
	if not line:
		break;

	parse = line.split()
	for index, xatoms in enumerate(L):
		atom1 = int(parse[2]); atom2 = int(parse[3])
		if atom1 in xatoms or atom2 in xatoms:
			E[index] += float(parse[5])

#angleinfo
while 1:
	line = angleinfo.readline()
	if not line:
		break;

	parse = line.split()
	for index, xatoms in enumerate(L):
		atom1 = int(parse[2]); atom2 = int(parse[3]); atom3 = int(parse[4])
		if atom1 in xatoms and atom2 in xatoms:
			E[index] += float(parse[6])
		elif atom2 in xatoms and atom3 in xatoms:
			E[index] += float(parse[6])
		elif atom1 in xatoms and atom3 in xatoms:
			E[index] += float(parse[6])

for index, value in enumerate(E):
	print(str(L[index][0]) + " " + "{0:8.5f}".format(value) + " " + str(xlnktype[L[index][0]]))
