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
from libposcar import modify_POSCAR

def pos_bombardment(pos, job, bomb_atoms, temp, nlevel, outfile):
    '''
    pos         POSCAR
    job         add     -> append atoms to POSCAR
                no add  -> select atoms in POSCAR
                
    bomb_atoms  a..  add atoms then they will be selected for bomb
                    aO20    add 20 O's
                s..  select atoms in POSCAR
                    sO2     select 2nd O group in atom list of POSCAR
                    sOn     select 1st O group and n atoms
    outfile     POSCAR.name
    '''
    
    #if re.match('s', bomb_atoms):
    modify_POSCAR(pos, job=job, matoms=bomb_atoms, option=temp, nlevel=nlevel, outf=outfile)
    #elif re.match('a', bomb_atoms):
    #    add_atoms(pos, bomb_atoms[1:])
    
    return 0

def main():
    parser = argparse.ArgumentParser(description="add atoms, vel block")
    parser.add_argument('poscar', help="poscar to be modified")
    parser.add_argument('-j', '--job', default='bomb', choices=['bomb','add','addbomb','zpe'], help="job of poscar changing")
    ### select existing atom or add atoms for bomb
    #parser.add_argument('-s', '--sel_atom', help="one atom species or index in POSCAR: Hf O1 Mo S O2 0 1 2 .. etc")
    parser.add_argument('-a', '--add_atoms',   help="atoms to be added")
    parser.add_argument('-ot', '--temp', type=float, default='298.15',  help="option of T for atom velocity")
    parser.add_argument('-l', '--nlevel', type=int, default=1,  help="atoms displaced in levels")
    gfname =  parser.add_mutually_exclusive_group()
    gfname.add_argument('-suf', '--suffix',     help="add suffix to outfile")
    gfname.add_argument('-o', '--outfile',      help='output POSCAR name')
    args = parser.parse_args()

    ### outfile name need to be passed
    if args.outfile:
        outfile = args.outfile
    elif args.suffix:
        outfile = args.poscar + args.suffix
    else:
        outfile = args.poscar + args.job

    #if 'bomb' in args.job:
    ### job = bomb or addbomb
    pos_bombardment(args.poscar, args.job, args.add_atoms, args.temp, args.nlevel, outfile)

    return 0

if __name__ == "__main__":
    main()
