#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys
import os
import bgf
import bgftools
import nutils as nu
import bgftools as bt
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

bottom = bt.atoms_average(b, 'atom.z', selection="'S_3a' in atom.ffType and atom.rNo == 1") # bottom
top = bt.atoms_average(b, 'atom.z', selection="'S_3b' in atom.ffType and atom.rNo == 2") # top

moved_atoms = ""

# round 3
for molecule in tqdm(mols):
    if len(molecule) != 3:
        continue
    cx, cy, cz = bgftools.getCom(b, ff, aNo_list=molecule)
    if abs(cz - top) < 0.5:
        for ano in molecule:
            atom = b.getAtom(ano)
            if 'I' in atom.chain:
                atom.z -= 2.0
                moved_atoms += "%d " % atom.aNo

# round 4
for molecule in tqdm(mols):
    if len(molecule) != 3:
        continue
    cx, cy, cz = bgftools.getCom(b, ff, aNo_list=molecule)
    if abs(cz - bottom) < 0.5:
        for ano in molecule:
            atom = b.getAtom(ano)
            if 'I' in atom.chain:
                atom.z += 2.0
                moved_atoms += "%d " % atom.aNo

# round 5
for molecule in tqdm(mols):
    if len(molecule) != 3:
        continue
    for ano in molecule:
        atom = b.getAtom(ano)
        if 'I' in atom.chain:
            if abs(atom.z - bottom) < 1.0:
                atom.z += 2.0
                moved_atoms += "%d " % atom.aNo

# round 6
for molecule in tqdm(mols):
    if len(molecule) != 3:
        continue
    for ano in molecule:
        atom = b.getAtom(ano)
        if 'I' in atom.chain:
            if abs(atom.z - top) < 1.0:
                atom.z -= 2.0
                moved_atoms += "%d " % atom.aNo


# round 1
for molecule in tqdm(mols):
    if len(molecule) != 3:
        continue
    cx, cy, cz = bgftools.getCom(b, ff, aNo_list=molecule)
    for ano in molecule:
        atom = b.getAtom(ano)
        if cz > top and atom.z < top: 
            atom.z += 1.0
            moved_atoms += "%d " % atom.aNo
        if cz < bottom and atom.z > bottom and "O" in atom.chain: 
            atom.z -= 2.0
            moved_atoms += "%d " % atom.aNo

# round 2
for molecule in tqdm(mols):
    if len(molecule) != 3:
        continue
    oz = 0
    for ano in molecule:
        atom = b.getAtom(ano)
        if "O" in atom.ffType:
            oz = atom.z
        else:
            if oz > bottom and atom.z < bottom:
                atom.z += 2.0
                moved_atoms += "%d " % atom.aNo
            elif oz > top and atom.z < top:
                atom.z += 2.0
                moved_atoms += "%d " % atom.aNo

print("atom %s moved." % moved_atoms)
b.REMARK.append("H atoms adjusted by GRA_avoid_bad_contacts.py")
b.saveBGF(sys.argv[2])
