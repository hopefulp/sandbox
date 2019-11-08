#!/home/joonho/anaconda3/bin/python

import argparse

from ase.io import read

def file_convert(fname, i_ext, o_ext):
    read(fname)


    return 0

def main():
    parser = argparse.ArgumentParser(description="file conversion ")
    parser.add_argument('fname', action='store_true', help="input file")
    parser.add_argument('-it','--input_ext', help="input file extension")
    parser.add_argument('-ot','--output_ext', help="output file extension")
    args = parser.parse_args()

    file_convert(args.fname, args.input_ext, args.output_ext)

if __name__ == "__main__":
    main()
