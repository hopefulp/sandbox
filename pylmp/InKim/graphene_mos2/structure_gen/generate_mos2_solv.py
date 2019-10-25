#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys
import os
import math
import bgf
import copy
import nutils as nu
import bgftools as bt
import scipy.spatial
import tqdm

usage = "%s bgf_file out_file ff_file layer_distance \n\tbgf_file: a mos2 single layer file" % sys.argv[0]
if len(sys.argv) < 5:
    print(usage)
    sys.exit(0)

bgf_file = sys.argv[1]
out_file = sys.argv[2]
ff_file = sys.argv[3]
layer_dist = float(sys.argv[4])

moslayer = bgf.BgfFile(bgf_file)

def get_average_z(mybgf, fftype):
    avg = 0.0; n = 0;
    for atom in mybgf.a:
        if fftype in atom.ffType:
            avg += atom.z
            n += 1
    return avg/float(n)

# residue name and residue number
for atom in moslayer.a:
    atom.rName = "MOS"
    atom.rNo = 1

# average z coord for Mo
avg_mo_z = bt.atoms_average(moslayer, 'atom.z', selection="'Mo' in atom.ffType")

# assign ffType
for atom in tqdm.tqdm(moslayer.a, ncols=120, desc='Assigning ffTypes'):
    if "S" in atom.ffType:
        if atom.z > avg_mo_z:
            atom.ffType = "S_3a"
        else:
            atom.ffType = "S_3b"
    
# remove infinite boundary
n_del = 0
for atom in tqdm.tqdm(moslayer.a, ncols=120, desc='Removing pbc bonds'):
    a = [atom.x, atom.y, atom.z]
    connected_ano = atom.CONECT
    for ano in connected_ano:
        atom2 = moslayer.getAtom(ano)
        a2 = [atom2.x, atom2.y, atom2.z]
        dist = nu.dist(a, a2)
        pbc_dist = nu.pbc_dist(a, a2, moslayer.CRYSTX[:3])
        if dist != pbc_dist:
            moslayer.disconnect(moslayer.a2i[atom.aNo], moslayer.a2i[atom2.aNo])
            n_del += 1
print("%d bonds are disconnected." % n_del)

# remove infinite boundary
n_del = 0
for atom in tqdm.tqdm(moslayer.a, ncols=120, desc='Removing pbc bonds 2'):
    a = [atom.x, atom.y, atom.z]
    connected_ano = atom.CONECT
    for ano in connected_ano:
        atom2 = moslayer.getAtom(ano)
        a2 = [atom2.x, atom2.y, atom2.z]
        dist = nu.dist(a, a2)
        pbc_dist = nu.pbc_dist(a, a2, moslayer.CRYSTX[:3])
        if dist != pbc_dist:
            moslayer.disconnect(moslayer.a2i[atom.aNo], moslayer.a2i[atom2.aNo])
            n_del += 1
print("%d bonds are disconnected." % n_del)

# calculate average z coord
avg_s3a_z = bt.atoms_average(moslayer, 'atom.z', selection="'S_3a' in atom.ffType")
avg_s3b_z = bt.atoms_average(moslayer, 'atom.z', selection="'S_3b' in atom.ffType")

# copy
moslayer2 = copy.deepcopy(moslayer)

# move
for atom in tqdm.tqdm(moslayer.a, ncols=120, desc='Translating'):
    atom.rName = "MOS"
    atom.rNo = 2
    atom.z += layer_dist + 2 * (avg_s3a_z - avg_mo_z)

# merge
result = moslayer.merge(moslayer2, True)
result.renumber()

'''
Overall structure is...

S_3a2 ----
Mo      :rNo 2
S_3b2 ----
        < layer_distance >
S_3a1 ----
Mo      :rNo 1
S_3b1 ----
'''

# save
result.FF = '/qcfs/noische/research/graphene/from_kjkwac/mos2/DREIDING2.21.ff_kkj20150330_mos2'
result.saveBGF("_temp.bgf")

# centerBGF
center_bgf_cmd = "~tpascal/scripts/centerBGF.pl -b %s -f %s -c com_center -s %s" % ("_temp.bgf", ff_file, "_temp.bgf")
nu.shutup(); os.system(center_bgf_cmd); nu.say()

# add water inside interlayer
result = bgf.BgfFile('_temp.bgf')

# add water to the system
addsolvent_cmd = "~tpascal/scripts/addSolvent.pl -i %s -f %s -n 'xy: +/- 20.0 z: +/- 20.0' -w spc -s %s" % ("_temp.bgf", ff_file, "_temp.bgf")
nu.shutup(); os.system(addsolvent_cmd); nu.say()

# empty center water
result = bgf.BgfFile('_temp.bgf')
min_mo_x1, max_mo_x1 = bt.atoms_minmax(result, 'atom.x', selection="'Mo' in atom.ffType and atom.rNo == 1")
min_mo_y1, max_mo_y1 = bt.atoms_minmax(result, 'atom.y', selection="'Mo' in atom.ffType and atom.rNo == 1")
avg_s3a_z1 = bt.atoms_average(result, 'atom.z', selection="'S_3a' in atom.ffType and atom.rNo == 1")
avg_s3b_z2 = bt.atoms_average(result, 'atom.z', selection="'S_3b' in atom.ffType and atom.rNo == 2")
all_molecules = bt.getMoleculeList(result)
del_list = []; mos2_coords = []
for atom in result.a:
    if "MOS" in atom.rName:
        mos2_coords.append([atom.x, atom.y, atom.z])

mos2_tree = scipy.spatial.KDTree(mos2_coords, leafsize=len(mos2_coords)+1)

for molecule in tqdm.tqdm(all_molecules, ncols=120, desc='Removing bad contacts'):
    if not len(molecule) == 3: continue # applies only to water

    cx, cy, cz = bt.getCom(result, ff_file=ff_file, aNo_list=molecule, silent=False)

    if min_mo_x1 < cx < max_mo_x1 and min_mo_y1 < cy < max_mo_y1 and avg_s3a_z1 < cz < avg_s3b_z2:
        for i in molecule: del_list.append(result.a2i[i])

    if mos2_tree.query_ball_point([cx, cy, cz], 2.0):
        for i in molecule: del_list.append(result.a2i[i])

result.delAtoms(del_list, silent=False)
result.renumber()
result.saveBGF("_temp.bgf")

# centerBGF
center_bgf_cmd = "~tpascal/scripts/centerBGF.pl -b %s -f %s -c com_center -s %s" % ("_temp.bgf", ff_file, out_file)
nu.shutup(); os.system(center_bgf_cmd); nu.say()

# renumber residue numbers
bt.renumberMolecules(out_file, out_file)

print("Done.")
