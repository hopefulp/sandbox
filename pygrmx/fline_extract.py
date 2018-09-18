#!/gpfs/home/joonho/anaconda3/bin/python
#/usr/bin/python

import argparse
import re
import os
import sys
from common import *

def cmd(fin, interval, fout, i_init, i_fin):

    with open(fin, 'r') as f:
        i=0
        for line in f:
            if i % interval == 0:
                print(line)

    return 0 


def main():
    parser = argparse.ArgumentParser(description='line extraction from file and plot')
    parser.add_argument( 'fname', help='input file')
    parser.add_argument( 'line_interval', type='int', help='line interval between extraction')
    parser.add_argument( '-i', '--init', default=0, type='int', help='start line number')
    parser.add_argument( '-f', '--final', default=1000000 type='int', help='end line number')
    parser.add_argument( '-o', '--outf', default='ext.dat', help='end line number')

    args = parser.parse_args()

    cmd(args.fname, args.line_interval, args.outf, args.init, args.final)
    return 0

if __name__ == "__main__":
    main()
