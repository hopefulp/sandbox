#!/home/noische/python

import sys
import os

import bgf
import bgftools

"""
Info)    Atoms: 5295
Info)    Bonds: 3858
Info)    Angles: 0  Dihedrals: 0  Impropers: 0  Cross-terms: 0
Info)    Bondtypes: 0  Angletypes: 0  Dihedraltypes: 0  Impropertypes: 0
Info)    Residues: 1630
Info)    Waters: 1629
Info)    Segments: 1
Info)    Fragments: 1630   Protein: 0   Nucleic: 0
"""

usage = """
BGF_info.py bgf_file ff_file
Prints information on BGF file. (mass, n_atoms, ff_types, ...)

"""
if len(sys.argv) < 3:
    print(usage)
    sys.exit(0)

bgf_file = sys.argv[1]
ff_file = sys.argv[2]

mybgf = bgf.BgfFile(bgf_file)

mw = bgftools.getMass(mybgf, ff_file)   # molecular weight
n_atoms = len(mybgf.a)  # n_atoms

d_fftypes = set()
for atom in mybgf.a:
    d_fftypes.add(atom.ffType)

print("Molecule mass: %f" % mw)
print("Number of molecules: %f" % n_atoms)
print("FF types in the molecule: %s" % str(list(d_fftypes)))
