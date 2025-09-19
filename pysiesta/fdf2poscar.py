#!/usr/bin/env python
from oorinano import *
from oorinano.calculator.siesta import readAtomicStructure  as read_geo
from oorinano.calculator.vasp   import writeAtomicStructure as write_geo
import sys, argparse
from common import f_ext, f_root


def fdf2poscar(fname, oformat, cell_size):
    '''
    input format:   fdf
    '''
    ext = f_ext(fname)
    froot = f_root(fname)

    atom = read_geo(fname)

    ### in case input is molecule, obtain cell info from argument
    #print(f"cell {atom.get_cell()}")
    if atom.get_cell() is None:
        print("Set cell for molecule")
        if type(cell_size) == float:
            va,vb,vc = [cell_size, 0, 0], [0, cell_size, 0], [0,0,cell_size]
        atom.set_cell([ va, vb, vc])


    write_geo(atom, f'{froot}.poscar')

def main():
    parser = argparse.ArgumentParser(
        description="Convert file format"
    )
    parser.add_argument("inf", help="input file")
    parser.add_argument('-of', '--out_format', default='POSCAR', help="output file format")
    parser.add_argument('-c', '--cell', default=15.0, type=float, help="cell size for molecule")

    args = parser.parse_args()

    fdf2poscar(args.inf, args.out_format, args.cell)


if __name__ == "__main__":
    main()

