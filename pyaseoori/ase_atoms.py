#!/home/joonho/anaconda3/bin/python

import argparse
from ase import Atoms, Atom, units
import ase.build as bd
from ase.build import bulk
from ase.visualize import view
from ase.io.formats import read, iread, write, string2index
from random import uniform
import numpy as np
import sys
import re
from ase.cluster import *
from common import list2str

def atoms_modify(str, atom_ind):
    atoms = io.read(str)


    for ind in atom_ind:
        atoms.arrays['positions'][ind] += 

def main():
    parser = argparse.ArgumentParser('ASE builder for structure')
    parser.add_argument('structure', help='atomic structure file')
    parser.add_argument('-j', '--job', help='translation, ')
    parser.add_argument('-al', '--atoms', nargs='*', help='atom indices from 0')
    parser.add_argument('-t', '--translate', default = [0.], nargs='*', help='x1.5, [1,0, 1.2], etc')
    parser.add_argument('-u', '--usage', action='store_true', help='shows usage and exit')
    args = parser.parse_args()

    if args.translate:
        


    atoms_modify(args.structure, args.atoms, translate)