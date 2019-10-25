#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys
import random
import bgf
import bgftools as bt
import nutils as nu

usage = """
%s bgf_file out_file selection n_del
""" % sys.argv[0]

if len(sys.argv) < 3:
    print(usage)
    sys.exit(0)

mybgf = bgf.BgfFile(sys.argv[1])
out_file = sys.argv[2]
selection = sys.argv[3]
n_del = int(sys.argv[4])

resid = set()
for atom in mybgf.a:
    if eval(selection):
        resid.add(atom.rNo)

resid = list(resid)
if len(resid) < n_del:
    nu.die("Not enough molecules to choose %d molecules to delete." % n_del)

rnos = random.sample(resid, n_del)

delist = []
for atom in mybgf.a:
    if atom.rNo in rnos:
        delist.append(mybgf.a2i[atom.aNo])

mybgf.delAtoms(delist)
mybgf.renumber()
mybgf.saveBGF(out_file)
bt.renumberMolecules(out_file, out_file)
