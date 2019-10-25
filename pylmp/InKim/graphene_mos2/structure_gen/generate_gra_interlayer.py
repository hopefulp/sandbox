#!/home/noische/python

import sys
import os
import math
import bgf
import copy
import nutils as nu
import bgftools as bt
import tqdm

usage = "%s bgf_file out_file ff_file layer_distance n_confined_water ext_water_margin" % sys.argv[0]
if len(sys.argv) < 7:
    print(usage)
    sys.exit(0)

bgf_file = sys.argv[1]
out_file = sys.argv[2]
ff_file = sys.argv[3]
layer_dist = float(sys.argv[4])
n_wat = int(sys.argv[5])
margin = float(sys.argv[6])

gralayer = bgf.BgfFile(bgf_file)

# residue name and residue number
for atom in gralayer.a:
    atom.rName = "GRA"
    atom.rNo = 1

# average z coord for GRA
avg_gra_z = bt.atoms_average(gralayer, 'atom.z', selection="'C_2G' in atom.ffType")

# copy
gralayer2 = copy.deepcopy(gralayer)

# move
for atom in tqdm.tqdm(gralayer2.a, ncols=120, desc='Translating'):
    atom.rName = "GRA"
    atom.rNo = 2
    atom.z += layer_dist

# merge
result = gralayer.merge(gralayer2, True)
result.renumber()

'''
Overall structure is...

C_2G    :rNo 2
        < layer_distance >
C_2G    :rNo 1
'''

# modify pbc
result.CRYSTX[2] = layer_dist + 2 * margin

# save
result.saveBGF("__temp.bgf")
#bt.renumberMolecules('__temp.bgf', '__temp.bgf')

# centerBGF
center_bgf_cmd = "~tpascal/scripts/centerBGF.pl -b %s -f '%s' -c com_center -s %s" % ("__temp.bgf", ff_file, "__temp.bgf")
nu.shutup(); os.system(center_bgf_cmd); nu.say()
#bt.renumberMolecules('__temp.bgf', '__temp.bgf')
mybgf = bgf.BgfFile('__temp.bgf')
avg_gra_z_bottom = bt.atoms_average(mybgf, 'atom.z', selection="'C_2G' in atom.ffType and atom.rNo == 1")
avg_gra_z_top = bt.atoms_average(mybgf, 'atom.z', selection="'C_2G' in atom.ffType and atom.rNo == 2")

# prepare water box
solvent = bgf.BgfFile('/home/noische/scripts/dat/WAT/spc_box.bgf')
solvent = bt.stress_cell(solvent, str(1/math.pow(2.0, 0.3333)))
solvent.saveBGF('_solvent.bgf')

addsolv_reg_cmd = "/qcfs/noische/scripts/BGF_addSolvent_region.py -b __temp.bgf -f '%s' -o __temp.bgf -m '0 0 %f' -M '%f %f %f' -n %d -t _solvent.bgf -v I" % (ff_file, avg_gra_z_bottom + 1.5, result.CRYSTX[0], result.CRYSTX[1], avg_gra_z_top - 1.5, n_wat)
print("** Adding confined water..")
print(addsolv_reg_cmd)
os.system(addsolv_reg_cmd)
bt.renumberMolecules('__temp.bgf', '__temp.bgf')

# add water outside graphene: bottom
addsolv_reg_cmd = "/qcfs/noische/scripts/BGF_addSolvent_region.py -b __temp.bgf -f '%s' -o __temp.bgf -M '%f %f %f' -v O" % (ff_file, result.CRYSTX[0], result.CRYSTX[1], avg_gra_z_bottom - 1.5)
print("** Adding bottom water..")
print(addsolv_reg_cmd)
os.system(addsolv_reg_cmd)
bt.renumberMolecules('__temp.bgf', '__temp.bgf')

# add water outside graphene: top
addsolv_reg_cmd = "/qcfs/noische/scripts/BGF_addSolvent_region.py -b __temp.bgf -f '%s' -m 'box box %f' -o %s -v O" % (ff_file, avg_gra_z_top + 1.5, out_file)
print("** Adding top water..")
print(addsolv_reg_cmd)
os.system(addsolv_reg_cmd)

# renumber residue numbers
bt.renumberMolecules(out_file, out_file)

print("Done.")
