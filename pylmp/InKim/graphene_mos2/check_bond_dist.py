import sys
import os
import numpy as np
import tqdm
import bgf
import nutils as nu
from pprint import *

usage = """Usage: %s bgf_file bond_delta

Shows a histogram of C_2G-C_2G bond length distribution.""" % os.path.basename(sys.argv[0])

if len(sys.argv) < 2:
    print(usage)
    sys.exit(0)

bgffile = bgf.BgfFile(sys.argv[1])
if sys.argv[2]:
    delta = float(sys.argv[2])
else:
    delta = 0.001

c = [atom for atom in bgffile.a if "C_2G" in atom.ffType]

dist = []
for atom in tqdm.tqdm(c, ncols=120, desc="Calculating"):
    x = [atom.x, atom.y, atom.z]
    for ano in atom.CONECT:
        atom2 = bgffile.getAtom(ano)
        y = [atom2.x, atom2.y, atom2.z]
        dist.append(nu.pbc_dist(x, y, bgffile.CRYSTX[:3]))

dist = np.array(dist)
if delta < (dist.max() - dist.min()):
    print(dist.min(), dist.max(), dist.max()-dist.min())
    bin = np.arange(dist.min(), dist.max(), delta)
    distr, _ = np.histogram(dist, bins=bin)
    norm_distr = distr / float(distr.sum())
    pprint(norm_distr)
else:
    print("average: %s, stdev: %s" % (dist.mean(), dist.std()))
