#!/usr/bin/env python3

import argparse
from ase.io import read, write
from ase.visualize import view

def ase_view(fname):
    # read the structure from XSF file
    atoms = read(fname, format='fdf')
    write('structure.cif', atoms)
    # open interactive ASE GUI
    view(atoms)

    return 0

def main():
    parser = argparse.ArgumentParser(
        description="read file format and view"
    )
    parser.add_argument( 'fname', help="input filename")
    parser.add_argument('-u', '--usage',   action='store_true', help='print usage')

    args = parser.parse_args()

    ase_view(args.fname)

    return 0

if __name__ == "__main__":
    main()