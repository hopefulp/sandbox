#!/home/joonho/anaconda3/bin/python

import argparse
#import sys
#import re
#import os
#import numpy as np

#from common import whereami
#import chem_space as cs
from pymatgen.core import Structure
from libposcar import parse_poscar


def pos_d22c(pos, job, outfile):
    '''
    pos         POSCAR
    job         d2c direct to cartesian
                c2d cartesian to direct
    outfile     POSCAR.name
    '''
    structure =  Structure.from_file(pos)

    line, cd = parse_poscar(pos, block='cd')
    pre_lines = parse_poscar(pos, block='pre')
    coord, natom = parse_poscar(pos, block='coord')

    if cd == 'D':
        cdline = 'Cartesian \n'
        if not outfile:
            outfile = pos + 'cart'
    elif cd == 'C':
        cdline = 'Direct \n'
        if not outfile:
            outfile = pos + 'drct'

    with open(outfile, 'w') as f:
        f.writelines(pre_lines)
        f.write(cdline)
        for i in range(structure.num_sites):
            f.writelines("{:14f} {:14f} {:14f}\n".format(structure.cart_coords[i][0],structure.cart_coords[i][1],structure.cart_coords[i][2]))  
    
    return 0

def main():
    parser = argparse.ArgumentParser(description="add atoms, vel block")
    parser.add_argument('poscar', help="poscar to be modified")
    parser.add_argument('-j', '--job', default='d2c', choices=['d2c'], help="direct to cartesian")
    gfname =  parser.add_mutually_exclusive_group()
    gfname.add_argument('-suf', '--suffix',     help="add suffix to outfile")
    gfname.add_argument('-o', '--outfile',      help='output POSCAR name')
    args = parser.parse_args()

    ### job = bomb or addbomb
    pos_d22c(args.poscar, args.job, args.outfile)

    return 0

if __name__ == "__main__":
    main()
