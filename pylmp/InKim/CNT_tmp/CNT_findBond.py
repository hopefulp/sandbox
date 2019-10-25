#!/home/noische/Enthought/Canopy_64bit/User/bin/python
import sys
import bgf, bgftools
import nutils as nu
import numpy as np
import scipy.spatial
from tqdm import tqdm

usage = """
CNT_findBond.py bgf_file out_file
"""

if len(sys.argv) < 3:
    print usage
    sys.exit(0)

mybgf = bgf.BgfFile(sys.argv[1])
pbc = mybgf.CRYSTX[:3]

l_cnt = []; 
l_watO = []; l_watH = []; 
l_mosS = []; l_mosMo = []; 
pairs = []

# Classification
for atom in mybgf.a:
    # graphitic
    if "C" in atom.ffType or "CNT" in atom.rName:
        l_cnt.append(atom)
    # Water
    elif "O" in atom.ffType or "OW" in atom.rName:
        l_watO.append(atom)
    elif "H" in atom.ffType or "HW" in atom.rName:
        l_watH.append(atom)
    # MoS2
    elif "Mo" in atom.ffType or "Mo" in atom.rName:
        l_mosMo.append(atom)
    elif "S" in atom.ffType or "S" in atom.rName:
        l_mosS.append(atom)
    else:
        pass;

print("len(l_cnt): %s" % len(l_cnt))
print("len(l_watO): %s" % len(l_watO))
print("len(l_mosMo): %s" % len(l_mosMo))

# CNT
if len(l_cnt):
    coords = []; anos = [];
    for atom in l_cnt:
        coords.append([atom.x, atom.y, atom.z])
        anos.append(atom.aNo)

    coords = np.array(coords)
    for atom in tqdm(l_cnt, ncols=80, miniters=100, desc="Analyzing C atoms"):
        #if len(atom.CONECT) == 3: continue;
        x = np.array([atom.x, atom.y, atom.z])
        tree = scipy.spatial.KDTree(coords, leafsize=len(coords)+1)
        dist, ndx = tree.query([x], k=4)  # the minimal distance is the atom itself
        index = np.where((dist >= 1.4) & (dist <= 1.5))
        for i in ndx[index]:
            #pairs.append([atom.aNo, anos[i]])
            mybgf.connectAtoms(atom.aNo, anos[i])

# water
if len(l_watO):
    coords = []; anos = [];
    for atom in l_watO:
        coords.append([atom.x, atom.y, atom.z])
        anos.append(atom.aNo)

    for atom in tqdm(l_watH, ncols=80, miniters=100, desc="Analyzing water  "):
        #if len(atom.CONECT) == 1: continue;
        x = np.array([atom.x, atom.y, atom.z])
        tree = scipy.spatial.KDTree(coords, leafsize=len(coords)+1)
        dist, ndx = tree.query([x], k=2)  # the minimal distance is the atom itself
        index = np.where((dist >= 0.9) & (dist <= 1.1))
        for i in ndx[index]:
            #pairs.append([atom.aNo, anos[i]])
            mybgf.connectAtoms(atom.aNo, anos[i])

'''
# MoS2
if len(l_mosS):
    coords = []; anos = [];
    for atom in l_mosS:
        coords.append([atom.x, atom.y, atom.z])
        anos.append(atom.aNo)

    for atom in tqdm(l_mosMo, ncols=80, miniters=100, desc="Analyzing MoS2   "):
        if len(atom.CONECT) == 6: continue;
        x = np.array([atom.x, atom.y, atom.z])
        tree = scipy.spatial.KDTree(coords, leafsize=len(coords)+1)
        dist, ndx = tree.query([x], k=2)  # the minimal distance is the atom itself
        index = np.where((dist >= 2.0) & (dist <= 2.6))
        for i in ndx[index]:
            #pairs.append([atom.aNo, anos[i]])
            mybgf.connectAtoms(atom.aNo, anos[i])
'''

# connect atoms
#for i in tqdm(pairs, ncols=80, miniters=100, desc="Connecting.....  "):
#    mybgf.connectAtoms(i[0], i[1])

# check
for atom in l_cnt:
    if len(atom.CONECT) != 3:
        nu.warn("Suspicious connectivity found on CNT C atom: " + str(atom.aNo))

for atom in l_watO:
    if len(atom.CONECT) != 2:
        nu.warn("Suspicious connectivity found on water O atom: " + str(atom.aNo))

for atom in l_mosMo:
    if len(atom.CONECT) != 6:
        nu.warn("Suspicious connectivity found on MoS2 Mo atom: " + str(atom.aNo))
mybgf.saveBGF(sys.argv[2])
