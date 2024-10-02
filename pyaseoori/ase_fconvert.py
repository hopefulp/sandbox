#!/home/joonho/anaconda3/bin/python

import argparse
from ase.io import read, write
import numpy as np
from common import f_ext, f_root
import sys, re

###        fname : format
formats = { 'POSCAR': 'vasp'}

def ase_convert(ifile, ofile, iff, off, movement, amount_m):
    ### Read
    try:
        atoms = read(ifile, format=iff)
    ### if fails read file format
    except:
        print(f"try input file format such as -if [vasp|vasp-out|vasp-xdatcar|vasp-xml|...]")
    ###### Output file format
    ### 1. extension is fformat
    if off:
        outff = off
    elif re.search('\.', ofile):
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

    try: 
        write(ofile, atoms, format=outff)
        print(f"{ofile} is written")
    except: 
        print(f"Error in writing {ofile}")

    return 0

def main():
    parser = argparse.ArgumentParser(description="read extxyz and write POSCAR  ")
    parser.add_argument('inf', help="input file")
    parser.add_argument('outf',  help="output file")
    parser.add_argument('-if','--inpff',  help="input file format")
    parser.add_argument('-of','--outff',  help="output file format")
    manipulate = parser.add_argument_group()
    manipulate.add_argument('-m', '--move', choices={'t','r'},  help="move molecule by translate|rotate")
    manipulate.add_argument('-ma', '--amount_move', default=0.5, help="amount of movement: translate w.r.t. cell size")
    args = parser.parse_args()

    ase_convert(args.inf, args.outf, args.inpff, args.outff, args.move, args.amount_move )

if __name__ == "__main__":
    main()
