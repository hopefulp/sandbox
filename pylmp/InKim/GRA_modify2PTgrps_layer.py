#!/home/noische/python

import sys
import os
import getopt
import shutil
import math

from numpy import pi
import numpy as np
import bgf
from grpfile import *
import nutils as nu

version = "160428"

def main(bgf_filename, grps_filename, group_no, l_bin):

    # initialization
    grps_original_filename = grps_filename + ".original"

    # keep original grps file at the first run
    if not os.path.exists(grps_original_filename):
        shutil.copy(grps_filename, grps_original_filename)

    # read atoms from bgf file
    mybgf = bgf.BgfFile(bgf_filename)
    pbc = mybgf.CRYSTX[:3]
    r_vdw_C = 2.45	# value from the gofr of trajectory
    r_vdw_C = 1.70	# value from the gofr of trajectory

    # graphene bottom
    print "%s: Calculating graphene z-coord.." % sys.argv[0]
    GRA1 = []; GRA2 = [];
    for atom in mybgf.a:
        if "GRA" in atom.rName and atom.rNo == 1:
            GRA1.append(atom.z)
        elif "GRA" in atom.rName and atom.rNo == 2:
            GRA2.append(atom.z)

    GRA1 = np.array(GRA1); z_GRA1 = np.mean(GRA1)   # graphene z-coords
    GRA2 = np.array(GRA2); z_GRA2 = np.mean(GRA2)   # graphene z-coords

    if l_bin != []:
        # reassign radius_bin
        l_bin[-1] = z_GRA2 - z_GRA1  # touch the last value
        r_atoms = dict()
        for index, i in enumerate(l_bin):
            z = z_GRA1 + i
            l_bin[index] = z
            r_atoms[z] = []

        count = 0
        # water position
        print "%s: Calculating water positions.. -- split by z-coord: %s" % (sys.argv[0], l_bin)
        WATER = [];
        for atom in mybgf.a:
            if "O" in atom.ffType and "WAT" in atom.rName and "I" in atom.chain:
                zcoord = np.array(atom.z)
                ind = np.digitize(zcoord, l_bin) # find a suitable bin
                r_atoms[l_bin[ind]].append(atom.aNo)    # add water mols to the bin
                r_atoms[l_bin[ind]].append(atom.CONECT[0])
                r_atoms[l_bin[ind]].append(atom.CONECT[1])

        for i in r_atoms.keys():
            r_atoms[i].sort()

        min_key = min(r_atoms.keys())
        r_atoms.pop(min_key, None)  # remove the smallest bin value

        # rewrite grps file
        print "%s: Creating grps file.." % sys.argv[0]
        print("\tGroup %d will be split into %d groups." % (group_no, len(l_bin)))
        g = grpfile(grps_original_filename)
        for key in r_atoms.keys():
            ub = l_bin.index(key)
            lb = ub - 1
            new_group_no = g.split_group(group_no, r_atoms[key])
            if new_group_no: print("A new group %d created from %d: range: %8.3f <= z <= %8.3f, natoms: %d" % (new_group_no, group_no, l_bin[lb], l_bin[ub], len(r_atoms[key])))

        #g.write(grps_filename, zip=True)
        g.write(grps_filename)
        
    print "%s: Done." % sys.argv[0]
    ### end of code


if __name__ == "__main__":

    usage = """
GRA_modify2PTgrps_layer.py -b bgf_file -g grps_file -G group -r layer_bin

- This script modifies Group Atoms in *grps file
- This script should be run within Tod's LAMMPS-2PT pbs script.
- layer_bin: e.g.) "2 8 10"
        the distance from the bottom layer (rName=="GRA" and rNo==1)
        the last value will be modified from the positions of carbon atoms in the BGF file.
"""

    if len(sys.argv) < 2:	
        print(usage);
        sys.exit(1)

    bgf_file = ""; grps_file = ""; layer_bin = ""; group = 0;

    options, args = getopt.getopt(sys.argv[1:], 'hb:g:r:G:', ['help', 'bgffile=', 'grpsfile=', 'bin=', 'group='])
    for option, value in options:
        if option in ('-h', '--help'):
            print(usage)
            sys.exit(0);
        elif option in ('-b', '--bgffile'):
            bgf_file = value
        elif option in ('-g', '--grpsfile'):
            grps_file = value
        elif option in ('-r', '--bin'):
            layer_bin = value
        elif option in ('-G', '--group'):
            group = int(value)
        

    layer_bin = [float(i) for i in layer_bin.split()]
    layer_bin.sort()

    main(bgf_file, grps_file, group, layer_bin)
