#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
import numpy as np
from common import dir_files

def extract_numeric_block(flines, i, j):
    return [list(map(float, line.strip().split())) for line in flines[i:j]]

def convert(job, inf, iline):
    ### Make 2D lines
    with open(inf, 'r') as f:
        flines = f.readlines()
    data_2d = extract_numeric_block(flines, iline[0], iline[1]+1)
    if job == 'D2to1':
        flattend = np.array(data_2d).flatten()
    print('\n'.join(map(str, flattend)))
    #print(f"length: {len(flattend)}")

    return 0 
    

def main():

    parser = argparse.ArgumentParser(description="display Usage for /mymplot  ")
    parser.add_argument('-j','--job', default='D2to1', help="data convert ")
    parser.add_argument('-i','--inf',  nargs='*', help="input file")
    parser.add_argument('-l','--iline', nargs=2,  type=int, help="line number in the file")
    parser.add_argument('-m','--math', default='diff', help="calculate difference ")
    args = parser.parse_args()

    convert(args.job, args.inf, args.iline, args.math)

if __name__ == "__main__":
    main()
