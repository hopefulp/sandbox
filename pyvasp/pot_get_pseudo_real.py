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

def get_pseudopot(infs, inp_type, iline, gmax, Lkspace, job, math, outfile):
    ### Make 2D lines

    D2_flist = []
    if outfile:
        fwrite = []
        outf = open(outfile, 'w')

    for i, inf in enumerate(infs):
        with open(inf, 'r') as f:
            flines = f.readlines()
        data_2d = extract_numeric_block(flines, iline[0], iline[1]+1)
        ### dimension of input potential
        if inp_type == 2:
            flattend = np.array(data_2d).flatten()
        nlen = len(list(flattend))
        
        if i == 0:
            D2_flist.append(list(flattend))
            nlen_old = nlen
            #print('\n'.join(map(str, flattend)))
        if i != 0 and nlen_old == len(list(flattend)):
            D2_flist.append(list(flattend))

    ncol2D = np.array(D2_flist).T
    N = ncol2D.shape[0]
    r_max = 10  # in Angstrom
    r = np.linspace( 0, r_max, N, endpoint=False)
    #delta_k = k[1] - k[0]
    #delta_x = 1/(N * delta_k)
    #x = np.arange(N) * delta_x
    if not Lkspace:
        for i, ncol in enumerate(ncol2D):
            print(f"{r[i]:>10.5f}", end=" ")
            s = f"{r[i]:>10.5f} "
            for col in ncol:
                print(f"{col:12.7e}", end=" ")
                s += f"{col:12.7e} "
            print()
            s += "\n"
            outf.write(s)
        print(f"{N} real-grid for pseudopotential to write to {outfile}")
    ### apply IFFT
    else :
        pass
        '''
        fft_y = np.fft.ifft(np.fft.ifftshift(ncol2D, axes=0), axis=0)
        y_real = np.real(fft_y)
        for i, ncol in enumerate(y_real):
            print(f"{x[i]:>10.5f}", end=" ")
            col_old = 0
            for col in ncol:
                print(f"{col:>15.7e}", end=" ")
                if col_old:
                    diff = col - col_old
                col_old = col
            print(f"{diff:>15.7e}")
        #print(f"length: {len(flattend)}")
        '''
    return 0 
    

def main():

    parser = argparse.ArgumentParser(description="Extract pseudopotential from POTCAR ")
    parser.add_argument('-i','--inf',  nargs='*', help="several POTCARs")
    parser.add_argument('-o','--outf',  default='pseudo.dat', help="pseudopotential outfile")
    parser.add_argument('-t','--in_type', default=2, type=int, help="convert 2D array to 1D")
    parser.add_argument('-l','--iline', nargs=2,  type=int, default=[57, 256], help="pseudopot part in POTCAR")
    parser.add_argument('-gm','--gmax', type=float, default=170.075338972103, help="to insert K-values in PSP at K-space")
    parser.add_argument('-m','--math', default='diff', help="calculate difference")
    parser.add_argument('-k','--kspace', action='store_true', help="convert k-space to R-space")
    parser.add_argument('-j','--job',  choices=['xscale','fft', 'ifft'], help="job will make 1D file w. fft, ifft, etc")
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

    get_pseudopot(args.inf, args.in_type, args.iline, args.gmax, args.kspace, args.job, args.math, args.outf)

if __name__ == "__main__":
    main()
