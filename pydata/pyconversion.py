#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os
import h5py

def conv_file(fname, ext):
    
    with h5py.File(fname, "r") as f:
        print("Top-level groups/datasets in file:")
        f.visititems(lambda name, obj: print(name, "->", type(obj)))
    return 0

def main():
    parser = argparse.ArgumentParser(description='To deal with file format')
    parser.add_argument('file', nargs='?', default='DOSCAR', help='read DOSCAR')
    parser.add_argument('-o', '--old', default='_o', help='save original DOSCAR to DOSCAR_o')
    args = parser.parse_args()
    
    conv_file(args.file, args.old) 

if __name__ == '__main__':
    main()
