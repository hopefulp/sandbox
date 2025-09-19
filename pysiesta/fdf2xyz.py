#!/usr/bin/env python
import argparse
from common import f_ext, f_root
import sys

oori_input_format = ['fdf']
output_format = ['xyz']

def convert_format(fname, oformat):
    '''
    input format:   fdf
    '''
    ext = f_ext(fname)
    froot = f_root(fname)
    if ext in oori_input_format:
        ### depending on input file format, import different library
        if ext == 'fdf':
            from oorinano.calculator.siesta import readAtomicStructure  as read_geo
            from oorinano import rw
            atoms = read_geo(fname)
        else:
            print("include more format in {oori_input_format}")
            sys.exit(1)

    if oformat in output_format:
        if oformat == 'xyz':
            rw.write_xyz(f'{froot}.xyz', atoms, comm=None, append=False)
    else:
        print("include more output file format in {output_format}")

def main():
    parser = argparse.ArgumentParser(
        description="Convert file format from fdf to any"
    )
    parser.add_argument("inf", help="input file")
    parser.add_argument('-of', '--out_format', default='xyz', help="output file format")

    args = parser.parse_args()

    convert_format(args.inf, args.out_format)


if __name__ == "__main__":
    main()
