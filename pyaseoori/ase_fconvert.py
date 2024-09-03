#!/home/joonho/anaconda3/bin/python

import argparse
from ase.io import read, write
import numpy as np
#from 

types = { 'POSCAR': 'vasp'}

def ase_convert(ifile, ofile, movement, amount_m):
    sys_bulk = read(ifile)
    ### check using view(import from ase.visualize) in ipython
    if movement:
        if movement=='t':
            length = sys_bulk.get_cell_lengths_and_angles()[0]
            arr = np.ones(sys_bulk.get_positions().shape)*length*amount_m
            sys_bulk.translate(arr)
        elif movement=='r':
            pass        # not ready


    #write(ofile, sys_bulk, format=types[ofile])
    write(ofile, sys_bulk, format='xsf')

    return 0

def main():
    parser = argparse.ArgumentParser(description="read extxyz and write POSCAR  ")
    parser.add_argument('-i','--in_file', help="input file")
    parser.add_argument('-o','--out_file',  help="output file")
    manipulate = parser.add_argument_group()
    manipulate.add_argument('-m', '--move', choices={'t','r'},  help="move molecule by translate|rotate")
    manipulate.add_argument('-ma', '--amount_move', default=0.5, help="amount of movement: translate w.r.t. cell size")
    args = parser.parse_args()

    ase_convert(args.in_file, args.out_file, args.move, args.amount_move )

if __name__ == "__main__":
    main()
