#!/home/joonho/anaconda3/bin/python

import argparse
from ase.io import read, write
import numpy as np
<<<<<<< HEAD
#from 
=======
from common import f_ext, f_root
import sys, re
>>>>>>> fb6711f6f651460c8fd75e3579a820fd55c5a3f2

###        fname : format
formats = { 'POSCAR': 'vasp'}

def ase_convert(ifile, ofile, off, movement, amount_m):
    ### Read
    atoms = read(ifile)

    ###### Output file format
    ### 1. extension is fformat
    if off:
        outff = off
    elif re.match('\.', ofile):
        outff = f_ext(ofile)     # output file format
    ### outfname has format in formats
    elif ofile in formats.keys():
        outff = formats[ofile]
    else:
        print(f"Can't recognize {ofile} file format")
        sys.exit(1)
    ### check using view(import from ase.visualize) in ipython
    if movement:
        if movement=='t':
            length = sys_bulk.get_cell_lengths_and_angles()[0]
            arr = np.ones(sys_bulk.get_positions().shape)*length*amount_m
            sys_bulk.translate(arr)
        elif movement=='r':
            pass        # not ready


<<<<<<< HEAD
    #write(ofile, sys_bulk, format=types[ofile])
    write(ofile, sys_bulk, format='xsf')
=======
    write(ofile, atoms, format=outff)
>>>>>>> fb6711f6f651460c8fd75e3579a820fd55c5a3f2

    return 0

def main():
    parser = argparse.ArgumentParser(description="read extxyz and write POSCAR  ")
    parser.add_argument('inf', help="input file")
    parser.add_argument('outf',  help="output file")
    parser.add_argument('-of','--outff',  help="output file format")
    manipulate = parser.add_argument_group()
    manipulate.add_argument('-m', '--move', choices={'t','r'},  help="move molecule by translate|rotate")
    manipulate.add_argument('-ma', '--amount_move', default=0.5, help="amount of movement: translate w.r.t. cell size")
    args = parser.parse_args()

    ase_convert(args.inf, args.outf, args.outff, args.move, args.amount_move )

if __name__ == "__main__":
    main()
