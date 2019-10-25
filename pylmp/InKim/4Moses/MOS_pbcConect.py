#!/home/noische/python
"""
MOS_pbcConect.py

Connects atoms near pbc in silica-amine system.
"""

import os, sys
import numpy as np
import bgf, bgftools
import nutils as nu

version = "150403"
usage = """
Usage: MOS_pbcConect.py amine_bgf (out_bgf)
"""

# script input
print(os.path.basename(sys.argv[0]) + " version " + str(version))
if len(sys.argv) < 2:
	print(usage); sys.exit(0)

bgf_file = sys.argv[1]
out_file = ""
if len(sys.argv) == 3:
	out_file = sys.argv[2]


# init
mybgf = bgf.BgfFile(bgf_file)
pbc = mybgf.CRYSTX[:3]


# loop over atoms 
# 1) Si-O pbc connection
pairs = [];
for atom in mybgf.a:
	if not "Si" in atom.ffType:
		continue;	# we need Si
	if len(atom.CONECT) == 4:
		continue;	# it should be undercoordinated

	x = np.array([atom.x, atom.y, atom.z])

	for atom2 in mybgf.a:
		if not "O" in atom2.ffType:
			continue;	# we need O
		if len(atom2.CONECT) == 2:
			continue;	# it should be undercoordinated

		y = np.array([atom2.x, atom2.y, atom2.z])
		d = bgftools.pbc_dist(x, y, pbc)

		# 1.62 ~ 1.63 A: from rdf of original BGF file
		if 1.6 < d < 1.7:
			pairs.append([atom.aNo, atom2.aNo])

print("    {0:<4d} atoms found for Si-O connection.".format(len(pairs)))


# 2) O-H pbc connection
pairs2 = [];
for atom in mybgf.a:
	if not "H" in atom.ffType:
		continue;

	x = np.array([atom.x, atom.y, atom.z])

	for atom2 in mybgf.a:
		if not "O" in atom2.ffType:
			continue;
		if len(atom2.CONECT) == 2:
			continue;

		y = np.array([atom2.x, atom2.y, atom2.z])
		d = bgftools.pbc_dist(x, y, pbc)

		if 1.62 < d < 1.63:
			# need to check distance when we use this later.. O-H bond length?
			pairs2.append([atom.aNo, atom2.aNo])

print("    {0:<4d} atoms found for O-H connection.".format(len(pairs2)))

pairs += pairs2	# merge


# connect atoms
for i in pairs:
	mybgf.connectAtoms(i[0], i[1])


# connectivity check-up
for atom in mybgf.a:
	if "Si" in atom.ffType:
		if len(atom.CONECT) != 4:
			nu.warn("Defect on Si connection found: " + str(atom.CONECTline()))
			x = np.array([atom.x, atom.y, atom.z])
			for ano in atom.CONECT:
				atom2 = mybgf.getAtom(ano)
				y = np.array([atom2.x, atom2.y, atom2.z])
				print bgftools.pbc_dist(x, y, pbc), atom2.ffType
	
	#print("    No connectivity error found")


# save
if not out_file:
	out_file = bgf_file[:-4] + ".pbc.bgf"

mybgf.saveBGF(out_file)
print("BGF saved to " + out_file)


print("Done.")
sys.exit(0)
