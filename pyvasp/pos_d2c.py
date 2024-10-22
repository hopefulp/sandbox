#!/home/joonho/anaconda3/bin/python
'''
add atom
add velocity block
'''
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
    line, cd = parse_poscar(pos, block='cd')
    pre_lines = parse_poscar(pos)

    if cd == 'D':
        cdline = 'Cartesian \n'
        if not outfile:
            outfile = pos + 'cart'
    elif cd == 'C':
        cdline = 'Direct \n'
        if not outfile:
            outfile = pos + 'drct'

    

    with open(pos, 'r') as fi and open(outfile, 'w') as fo:
        lines = fi.readlines()
        for i, line in enumerate(lines):
            if i < 8:
                fo.write(line)
            if i == 8:
                if 
    #elif re.match('a', bomb_atoms):
    #    add_atoms(pos, bomb_atoms[1:])
    
    return 0

def main():
    parser = argparse.ArgumentParser(description="add atoms, vel block")
    parser.add_argument('poscar', help="poscar to be modified")
    parser.add_argument('-j', '--job', default='d2c', choices=['d2c'], help="direct to cartesian")
    gfname =  parser.add_mutually_exclusive_group()
    gfname.add_argument('-suf', '--suffix',     help="add suffix to outfile")
    gfname.add_argument('-o', '--outfile',      help='output POSCAR name')
    args = parser.parse_args()

    ### outfile name need to be passed
    if args.outfile:
        outfile = args.outfile
    elif args.suffix:
        outfile = args.poscar + args.suffix


    #if 'bomb' in args.job:
    ### job = bomb or addbomb
    pos_d22c(args.poscar, args.job, outfile)

    return 0

if __name__ == "__main__":
    main()
