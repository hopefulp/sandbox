#!/home/noische/python
import sys
import bgf, bgftools
import nutils as nu
import numpy as np

usage = """
CNT_findBond.py bgf_file out_file
"""

if len(sys.argv) < 3:
    print usage
    sys.exit(0)

mybgf = bgf.BgfFile(sys.argv[1])
pbc = mybgf.CRYSTX[:3]

l_cnt = []; l_watO = []; l_watH = []; l_mosMo = []; l_mosS = []
pairs = []
print "sweeping"
for atom in mybgf.a:
    if "C" in atom.ffType or "CNT" in atom.rName:
        l_cnt.append(atom)
    elif "O" in atom.ffType or "OW" in atom.rName:
        l_watO.append(atom)
    elif "H" in atom.ffType or "HW" in atom.rName:
        l_watH.append(atom)
    elif "S" in atom.ffType:
        l_mosS.append(atom)
    elif "Mo" in atom.ffType:
        l_mosMo.append(atom)
    else:
        pass;

'''
# CNT
print "Reading CNT"
n_prog_total = len(l_cnt)
n_prog = 0
for atom1 in l_cnt:
    n_prog += 1
    x = np.array([atom1.x, atom1.y, atom1.z])
    sys.stdout.write("\r" + "Progress " + str(n_prog) + " / " + str(n_prog_total))
    sys.stdout.flush()
    for atom2 in l_cnt:
        y = np.array([atom2.x, atom2.y, atom2.z])
        d = bgftools.pbc_dist(x, y, pbc)
        dx = abs(atom.x - atom2.x)
        dy = abs(atom.y - atom2.y)
        dz = abs(atom.z - atom2.z)

        if 1.4 < d < 1.5:
            pairs.append([atom1.aNo, atom2.aNo])

# water
print "Reading WAT"
n_prog_total = len(l_watO)
n_prog = 0
for atom1 in l_watO:
    n_prog += 1
    x = np.array([atom1.x, atom1.y, atom1.z])
    sys.stdout.write("\r" + "Progress " + str(n_prog) + " / " + str(n_prog_total))
    sys.stdout.flush()
    for atom2 in l_watH:
        y = np.array([atom2.x, atom2.y, atom2.z])
        d = bgftools.pbc_dist(x, y, pbc)

        if 0.9 < d < 1.1:
            pairs.append([atom1.aNo, atom2.aNo])
'''

# mos2
print "Reading MOS"
n_prog_total = len(l_mosMo)
n_prog = 0
for atom1 in l_mosMo:
    n_prog += 1
    x = np.array([atom1.x, atom1.y, atom1.z])
    sys.stdout.write("\r" + "Progress " + str(n_prog) + " / " + str(n_prog_total))
    sys.stdout.flush()
    for atom2 in l_mosS:
        y = np.array([atom2.x, atom2.y, atom2.z])
        d = bgftools.pbc_dist(x, y, pbc)

        if 2.3 < d < 2.5:
            #pairs.append([atom1.aNo, atom2.aNo])
            mybgf.connectAtoms(atom1.aNo, atom2.aNo)

# connect atoms
print "connecting"
for i in pairs:
    mybgf.connectAtoms(i[0], i[1])

# check
'''
for atom in l_cnt:
    if len(atom.CONECT) != 3:
        nu.warn("Suspicious connectivity found on CNT C atom: " + str(atom.aNo))

for atom in l_watO:
    if len(atom.CONECT) != 2:
        nu.warn("Suspicious connectivity found on water O atom: " + str(atom.aNo))
'''

for atom in l_mosMo:
    if len(atom.CONECT) != 6:
        nu.warn("Suspicious connectivity found on CNT C atom: " + str(atom.aNo))

mybgf.saveBGF(sys.argv[2])
