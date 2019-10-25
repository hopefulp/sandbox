#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys
import os
import bgf
import bgftools
import nutils as nu
from tqdm import tqdm

usage = """
avoid.py bgffile outfile fffile
"""
if len(sys.argv) < 2:
    print usage
    sys.exit(0)

print(sys.argv)

b = bgf.BgfFile(sys.argv[1])
if len(sys.argv) <= 3:
    ff = os.environ['FF_CNT'].replace("'", "")
else:
    ff = sys.argv[3]

mols = bgftools.getMoleculeList(b)

# GRA 1 : bottom sheet
# GRA 2 : top sheet
top = 0.0; bottom = 0.0;
for atom in b.a:
    if "GRA" in atom.rName:
        if atom.rNo == 1:
            bottom = atom.z
        elif atom.rNo == 2:
            top = atom.z

if top and bottom:
    print("Found two graphene sheets.")
else:
    nu.die("Two graphene sheets with resname 1 (bottom) or 2 (top) found. Exiting.")

# round 1
moved_atoms = ""
for molecule in tqdm(mols):
    if len(molecule) != 3:
        continue
    cx, cy, cz = bgftools.getCom(b, ff, aNo_list=molecule)
    for ano in molecule:
        atom = b.getAtom(ano)
        if cz > top and atom.z < top: 
            atom.z += 1
            moved_atoms += "%d " % atom.aNo
        if cz < bottom and atom.z > bottom and "O" in atom.chain: 
            atom.z -= 2
            moved_atoms += "%d " % atom.aNo

# round 2
for molecule in tqdm(mols):
    oz = 0
    for ano in molecule:
        atom = b.getAtom(ano)
        if "O" in atom.ffType:
            oz = atom.z
        else:
            if oz > bottom and atom.z < bottom:
                atom.z += 1.0
                moved_atoms += "%d " % atom.aNo
            elif oz > top and atom.z < top:
                atom.z += 1.0
                moved_atoms += "%d " % atom.aNo

print("atom %s moved." % moved_atoms)
b.REMARK.append("H atoms adjusted by GRA_avoid_bad_contacts.py")
b.saveBGF(sys.argv[2])
