#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os

def modify_dos(fname, ext):
    return x

def main():
    parser = argparse.ArgumentParser(description='To get rid of abnormal DOS at start energy')
    parser.add_argument('file', nargs='?', default='DOSCAR', help='read DOSCAR')
    parser.add_argument('-o', '--old', default='_o', help='save original DOSCAR to DOSCAR_o')
    args = parser.parse_args()
    
    modify_doscar(args.file, args.old) 

if __name__ == '__main__':
    main()
