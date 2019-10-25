#!/home/noische/python

import sys
import os
import time
import shutil
import math

from numpy import pi
import numpy as np
import bgf
import CNT
from LAMMPS_trj2bgf import *
import nutils as nu

version = "160303"

def main(bgf_filename, grps_filename, r_bin):

    # initialization
    grps_original_filename = grps_filename + ".original"

    if r_bin == "":
        radius_bin = []
    else:
        radius_bin = r_bin

    n_group = len(radius_bin) + 2  # CNT, total, 0~8, 8~15, 15~22

    # keep original grps file at the first run
    if not os.path.exists(grps_original_filename):
        shutil.copy(grps_filename, grps_original_filename)

    # read atoms from bgf file
    mybgf = bgf.BgfFile(bgf_filename)
    pbc = mybgf.CRYSTX[:3]
    r_vdw_C = 2.45	# value from the gofr of trajectory

    # CNT center
    print "%s: Calculating CNT center.." % sys.argv[0]
    CNT = []; cnt_atoms = []; water_atoms = [];
    for atom in mybgf.a:
        if "NT" in atom.rName:
            CNT.append([atom.x, atom.y])
            cnt_atoms.append(atom.aNo)
        elif "WAT" in atom.rName or "OW" in atom.ffType or "HW" in atom.ffType:
            water_atoms.append(atom.aNo)
    cnt_atoms.sort()
    water_atoms.sort()

    CNT = np.array(CNT); n_CNT = len(CNT); n_water = len(water_atoms) / 3
    x_CNT, y_CNT = np.mean(CNT, axis=0)
    CNT_center = np.array([x_CNT, y_CNT])   # CNT center
    CNT_height = mybgf.CRYSTX[2]        # CNT height from pbc
    CNT_radius = np.mean(np.sqrt(((CNT-CNT_center)**2).sum(axis=-1)))

    if radius_bin != []:
        # reassign radius_bin
        radius_bin[-1] = CNT_radius - r_vdw_C  # touch the last value
        r_atoms = dict()
        for i in radius_bin:
            r_atoms[i] = []

        # water distance
        print "%s: Calculating water distance.." % sys.argv[0]
        WATER = [];
        for atom in mybgf.a:
            if "O" in atom.ffType:
                coord = np.array([atom.x, atom.y])
                dist = np.sqrt(((coord - CNT_center)**2).sum(axis=-1))  # O distance from CNT center
                ind = np.digitize(dist, radius_bin) # find a suitable bin
                r_atoms[radius_bin[ind]].append(atom.aNo)    # add water mols to the bin
                r_atoms[radius_bin[ind]].append(atom.CONECT[0])
                r_atoms[radius_bin[ind]].append(atom.CONECT[1])

        for i in r_atoms.keys():
            r_atoms[i].sort()

    # effective volume
    total_volume = pbc[0] * pbc[1] * pbc[2]
    vacuum_volume = total_volume - ( (CNT_radius + r_vdw_C)**2 * math.pi * pbc[2] )
    cnt_volume = ((CNT_radius + r_vdw_C)**2 - (CNT_radius - r_vdw_C)**2) * math.pi * pbc[2]
    inner_volume = 0.0; group_volume = []; r = [0] + radius_bin # 0 8 15 x
    for i in range(len(radius_bin)):
        vol = (r[i+1]**2 - r[i]**2) * math.pi * pbc[2]
        group_volume.append(vol)
        inner_volume += vol
    
    if radius_bin == []:
        inner_volume = (CNT_radius - r_vdw_C)**2 * math.pi * pbc[2]

    cnt_volume = total_volume - inner_volume    #

    # write grps file
    print "%s: Creating grps file.." % sys.argv[0]
    line = ""
    line += "Total Groups: " + str(n_group) + "\n"
    # Group Atoms: CNT group
    line += "Group %d Atoms %d" % (1, n_CNT) + "\n"
    line += ' '.join((('%i - %i' % r) if len(r) == 2 else '%i' % r) for r in nu.range_extract(cnt_atoms)) + "\n"
    # Group Atoms: total water group
    line += "Group %d Atoms %d" % (2, len(water_atoms)) + "\n"
    line += ' '.join((('%i - %i' % r) if len(r) == 2 else '%i' % r) for r in nu.range_extract(water_atoms)) + "\n"
    # Group Atoms: radial water groups
    for index, i in enumerate(radius_bin):
        grpid = index + 3   # 1 + cnt + total_water
        line += "Group %d Atoms %d" % (grpid, len(r_atoms[i])) + "\n"
        ##line += ' '.join((('%i - %i' % r) if len(r) == 2 else '%i' % r) for r in nu.range_extract(r_atoms[i])) + "\n" # failed
        line += ' '.join("%d " % i for i in r_atoms[i]) + "\n"

    # Constraints
    line += "Constraints" + "\n"
    line += "0 %d " % (len(water_atoms))
    for index, i in enumerate(radius_bin):
        line += "%d " % (len(r_atoms[i]))
    line += "\n"

    # RotationalSymmetryNumber
    line += "RotationalSymmetryNumber" + "\n"
    line += "1 2 "
    for index, i in enumerate(radius_bin):
        line += "2 "
    line += "\n"

    # LinearMoleculeFlag
    line += "LinearMoleculeFlag" + "\n"
    line += "0 0 "
    for index, i in enumerate(radius_bin):
        line += "0 "
    line += "\n"

    # GroupVolume
    line += "GroupVolume" + "\n"
    line += "%15.6f %15.6f " % (cnt_volume, inner_volume)
    line += ''.join("%15.6f " % i for i in group_volume) + "\n"

    gf = open(grps_filename, 'w')	# reopen
    gf.writelines(line)
    gf.close()

    print "%s: Done." % sys.argv[0]
    ### end of code


if __name__ == "__main__":

    usage = """
CNT_modify2PTgrps_radius.py bgf_file grps_file radial_bin

- This script modifies Group Atoms in *grps file
- Originally intended to use for (36, 36) CNT with four radial sections. (now works for various radial bins)
- This script should be run within Tod's LAMMPS-2PT pbs script.
- radial_bin: "2 8 10" -> the last value will be modified from the positions of carbon atoms in the BGF file.
"""

    if len(sys.argv) < 2:	
        print(usage);
        sys.exit(1)

    bgf_file = ""; grps_file = ""; radial_bin = "";

    options, args = getopt.getopt(sys.argv[1:], 'hb:g:r:', ['help', 'bgffile=', 'grpsfile=', 'bin='])
    for option, value in options:
        if option in ('-h', '--help'):
            print(usage)
            sys.exit(0);
        elif option in ('-b', '--bgffile'):
            bgf_file = value
        elif option in ('-g', '--grpsfile'):
            grps_file = value
        elif option in ('-r', '--bin'):
            radial_bin = value

    radial_bin = [int(i) for i in list(radial_bin) if i != " "]
    radial_bin.sort()

    main(bgf_file, grps_file, radial_bin)
