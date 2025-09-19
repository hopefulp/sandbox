#!/usr/bin/env python
import argparse
from common import f_ext, f_root
import sys

oori_input_format = ['fdf']
ase_input_format = ['xyz','cif','POSCAR','xsf']
output_format = ['xyz']

def convert_format(fname, oformat, Lrevert):
    '''
    input format:   fdf
    '''
    if '.' in fname:
        ext     = f_ext(fname)
        froot   = f_root(fname)
    else:
        print(f"input should have extension {fname}")
        sys.exit(11)

    if Lrevert:
        tmp     = ext
        froot   = ext
        ext     = tmp
    outf    = froot + f'.{oformat}'
    if ext in oori_input_format:
        ### depending on input file format, import different library
        from oorinano.calculator.siesta import readAtomicStructure  as read
        from oorinano import rw as write
    elif ext in ase_input_format:
        from ase.io import read, write
    else:
        print("include more format in {oori_input_format} or {ase_input_format}")
        sys.exit(1)
    atoms = read(fname)

    ### write
    if oformat in output_format:
        if ext in oori_input_format:
            if oformat == 'xyz':
                write.write_xyz(f'{froot}.xyz', atoms, comm=None, append=False)
    elif ext in ase_input_format:
        write(outf, atoms)
    else:
        print("include more output file format in {output_format}")

def main():
    parser = argparse.ArgumentParser(
        description="Generate bond dissociation geometries for diatomic molecule."
    )
    parser.add_argument("inf", help="input file")
    parser.add_argument('-of', '--out_format', default='xyz', help="output file format")
    parser.add_argument('-r', '--revert', action='store_true', help="extension is first e.g. POSCAR.name")

    args = parser.parse_args()

    convert_format(args.inf, args.out_format, args.revert)


if __name__ == "__main__":
    main()
