#!/opt/applic/epd/bin/python

import bgf
import sys

bgf1 = bgf.BgfFile(sys.argv[1])

for atom in bgf1.a:
	if "C_" in atom.ffType:
		n_H = [];
		for ano in atom.CONECT:
			if bgf1.getAtom(ano).is_hydrogen(): n_H.append(ano)
		print('atom No: %s  n_H: %s' % (atom.aNo, n_H))
		if len(n_H) != 3:
			bgf1.delAtoms(n_H)
print(bgf1)
