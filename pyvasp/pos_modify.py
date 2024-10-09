#!/home/joonho/anaconda3/bin/python
'''
add atom
add velocity block
'''
import argparse
import sys
import re
import os
import numpy as np

from common import whereami
import chem_space as cs
from libposcar import modify_POSCAR

def pos_addatoms(pos, add_atoms):
    '''
    add add_atoms at the end of POSCAR
    and those are bombardment
    '''

    



    return 0



def pos_bombardment(pos, bomb_atoms, temp, outfile):
    '''
    pos         POSCAR
    bomb_atoms  a..  add atoms then they will be selected for bomb
                    aO20    add 20 O's
                s..  select atoms in POSCAR
                    sO2     select 2nd O group in atom list of POSCAR
                    sOn     select 1st O group and n atoms
    outfile     POSCAR.name
    '''
    
    if re.match('s', bomb_atoms):
        modify_POSCAR(pos, job='vel', sel_atom=bomb_atoms[1:], option=temp, outf=outfile)
    elif re.match('a', bomb_atoms):
        add_atoms(pos, bomb_atoms[1:])
    
    return 0

def main():
    parser = argparse.ArgumentParser(description="add atoms, vel block")
    parser.add_argument('poscar', help="poscar to be modified")
    parser.add_argument('-j', '--job', default='bomb', choices=['bomb','add','addbomb'], help="job of poscar changing")
    ### select existing atom or add atoms for bomb
    gatoms = parser.add_mutually_exclusive_group()
    gatoms.add_argument('-s', '--sel_atom', help="one atom species or index in POSCAR: Hf O1 Mo S O2 0 1 2 .. etc")
    gatoms.add_argument('-a', '--add_atoms',   help="atoms to be added")
    parser.add_argument('-ot', '--temp', type=float, default='298.15',  help="option of T for atom velocity")
    gfname =  parser.add_mutually_exclusive_group()
    gfname.add_argument('-suf', '--suffix',     help="atoms to be added")
    gfname.add_argument('-o', '--outfile',      help='output POSCAR name')
    args = parser.parse_args()

    ### outfile name need to be passed
    if args.outfile:
        outfile = args.outfile
    elif args.suffix:
        outfile = args.poscar + args.suffix
    else:
        outfile = args.poscar + args.job

    if args.sel_atom:
        atoms = 's' + args.sel_atom
    elif args.add_atoms:
        atoms = 'a' + args.add_atoms

    if 'bomb' in args.job:
        ### job = bomb or addbomb
        pos_bombardment(args.poscar, atoms, args.temp, outfile)
    '''
    ### 2D scattered atom coordinates are required
    elif args.job == 'add':
        pos_add_atom(
    '''
    return 0

if __name__ == "__main__":
    main()
