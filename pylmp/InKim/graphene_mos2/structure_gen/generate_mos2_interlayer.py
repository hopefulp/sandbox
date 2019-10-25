#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys
import os
import math
import bgf
import copy
import nutils as nu
import bgftools as bt
import tqdm

usage = "%s bgf_file out_file ff_file layer_distance n_confined_water ext_water_margin (no if no_ext_water)" % sys.argv[0]
if len(sys.argv) < 7:
    print(usage)
    sys.exit(0)

bgf_file = sys.argv[1]
out_file = sys.argv[2]
ff_file = sys.argv[3]
layer_dist = float(sys.argv[4])
n_wat = int(sys.argv[5])
margin = float(sys.argv[6])

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
    
# calculate average z coord
avg_s3a_z = bt.atoms_average(moslayer, 'atom.z', selection="'S_3a' in atom.ffType")
avg_s3b_z = bt.atoms_average(moslayer, 'atom.z', selection="'S_3b' in atom.ffType")

# copy
moslayer2 = copy.deepcopy(moslayer)

# move
for atom in tqdm.tqdm(moslayer2.a, ncols=120, desc='Translating'):
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

# modify pbc
result.CRYSTX[2] = layer_dist + 4 * (avg_s3a_z - avg_mo_z) + 2 * margin

# save
result.saveBGF("_temp.bgf")
#bt.renumberMolecules('_temp.bgf', '_temp.bgf')

# centerBGF
center_bgf_cmd = "~tpascal/scripts/centerBGF.pl -b %s -f %s -c com_center -s %s" % ("_temp.bgf", ff_file, "_temp.bgf")
nu.shutup(); os.system(center_bgf_cmd); nu.say()
#bt.renumberMolecules('_temp.bgf', '_temp.bgf')
mybgf = bgf.BgfFile('_temp.bgf')

# add water inside interlayer
avg_s3a_z1 = bt.atoms_average(mybgf, 'atom.z', selection="'S_3a' in atom.ffType and atom.rNo == 1")
avg_s3b_z1 = bt.atoms_average(mybgf, 'atom.z', selection="'S_3b' in atom.ffType and atom.rNo == 1")
avg_s3a_z2 = bt.atoms_average(mybgf, 'atom.z', selection="'S_3a' in atom.ffType and atom.rNo == 2")
avg_s3b_z2 = bt.atoms_average(mybgf, 'atom.z', selection="'S_3b' in atom.ffType and atom.rNo == 2")

# prepare water box
solvent = bgf.BgfFile('/home/noische/scripts/dat/WAT/spc_box.bgf')
solvent = bt.stress_cell(solvent, str(1/math.pow(1.7, 0.3333)))
solvent.saveBGF('_solvent.bgf')

addsolv_reg_cmd = "/qcfs/noische/scripts/BGF_addSolvent_region.py -b _temp.bgf -f %s -o _temp.bgf -m '0 0 %f' -M '%f %f %f' -n %d -t _solvent.bgf -v I" % (ff_file, avg_s3a_z1 + 1.5, result.CRYSTX[0], result.CRYSTX[1], avg_s3b_z2 - 1.0, n_wat)
print("** Adding confined water..")
print(addsolv_reg_cmd)
os.system(addsolv_reg_cmd)
bt.renumberMolecules('_temp.bgf', '_temp.bgf')

if len(sys.argv) == 8:
    os.system("mv _temp.bgf %s" % out_file)
    sys.exit(0)

# add water outside mos2: bottom
addsolv_reg_cmd = "/qcfs/noische/scripts/BGF_addSolvent_region.py -b _temp.bgf -f %s -o _temp.bgf -M '%f %f %f' -v O" % (ff_file, result.CRYSTX[0], result.CRYSTX[1], avg_s3b_z1 - 0.5)
print("** Adding bottom water..")
print(addsolv_reg_cmd)
os.system(addsolv_reg_cmd)
bt.renumberMolecules('_temp.bgf', '_temp.bgf')

# add water outside mos2: top
addsolv_reg_cmd = "/qcfs/noische/scripts/BGF_addSolvent_region.py -b _temp.bgf -f %s -m '0 0 %f' -o %s -v O" % (ff_file, avg_s3a_z2 + 1.5, out_file)
print("** Adding top water..")
print(addsolv_reg_cmd)
os.system(addsolv_reg_cmd)

# renumber residue numbers
bt.renumberMolecules(out_file, out_file)

print("Done.")
