#!/home/joonho/anaconda3/bin/python

import argparse
import re
from ase.io import read, write

def ase_convert(inf, outtype, job, outfname):
    '''
    Using ASE, read POSCAR/CONTCAR write POSCAR w. c/d
    inf         input filename
    outtype     output file type
    job         convert (default)
    outfile     outfile (default)
    '''
    ### multi frames
    Lmulti = 0
    if inf == 'OUTCAR':
        Lmulti = 1
    if Lmulti:
        atoms = read(inf, index=":")
    ### single frames
    else:
        atoms = read(inf)

    write(outfname, atoms)

    print(f"{inf} is written to {outfname}")
    return 0

def main():
    parser = argparse.ArgumentParser(description="add atoms, vel block")
    parser.add_argument( 'inf', help="file with vasp type")
    parser.add_argument('-t', '--type', help="output file type")
    parser.add_argument('-j', '--job', default='convert', choices=['convert'], help="convert file type")
    gfname =  parser.add_mutually_exclusive_group()
    gfname.add_argument('-suf', '--suffix',     help="add suffix to outfile")
    gfname.add_argument('-o', '--outfile', default='outfile',  help='output POSCAR name')
    args = parser.parse_args()

    if args.suffix:
        outf = args.inf + args.suffix
    else:
        outf = args.outfile

    ### job = bomb or addbomb
    ase_convert(args.inf,  args.type, args.job, outf)

    return 0

if __name__ == "__main__":
    main()
