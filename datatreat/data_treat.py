#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
import sys
import numpy as np
from common import dir_files
from my_data import extract_numeric_block
#def extract_numeric_block(flines, i, j):
#    return [list(map(float, line.strip().split())) for line in flines[i:j]]

def convert(job, infs, iline, math):
    ### Make 2D lines

    D2_flist = []

    for i, inf in enumerate(infs):
        with open(inf, 'r') as f:
            flines = f.readlines()
        data_2d = extract_numeric_block(flines, iline[0], iline[1]+1)
        if job == 'D2to1':
            flattend = np.array(data_2d).flatten()
        nlen = len(list(flattend))
        
        if i == 0:
            D2_flist.append(list(flattend))
            nlen_old = nlen
            #print('\n'.join(map(str, flattend)))
        if i != 0 and nlen_old == len(list(flattend)):
            D2_flist.append(list(flattend))

    if len(infs) == 2:
        ncol2D = np.array(D2_flist).T
        for i, ncol in enumerate(ncol2D):
            #print(f"{i:>10d} {ncol[0]:12.7e} {ncol[1]:12.7e} {ncol[1]-ncol[0]:12.7e}")
            print(f"{ncol[0]:12.7e} {ncol[1]:12.7e}")
            
            
        #print(f"length: {len(flattend)}")

    return 0 
    

def main():

    parser = argparse.ArgumentParser(description="display Usage for /mymplot  ")
    parser.add_argument('-j','--job', default='D2to1', help="data convert ")
    parser.add_argument('-i','--inf',  nargs='*', help="input file")
    parser.add_argument('-l','--iline', nargs=2,  type=int, help="line number in the file")
    parser.add_argument('-m','--math', default='diff', help="calculate difference ")
    parser.add_argument('-u','--usage', action='store_true', help="show usage ")
    args = parser.parse_args()

    if args.usage:
        print("\
            \n\t\tdata_treat.py -i POTCAR2.H/0/POTCAR POTCAR2.H/0.1/POTCAR -l 57 256 > diff1.txt\
            \n\t\t(plot)\
            \n\t\t    gnu_1f_nc.sh diff1.txt 2 3 4\
            \n\t\t    gnu_1f_nc.sh diff1.txt 4\
            \n\t\t\
            \n")
        sys.exit(0)

    convert(args.job, args.inf, args.iline, args.math)

if __name__ == "__main__":
    main()
