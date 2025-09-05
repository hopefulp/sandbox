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

def get_pseudopot(infs, inp_type, iline, gmax, job, math, outfile):
    ### Make 2D lines

    if outfile:
        fwrite = []
        outf = open(outfile, 'w')

    v_g_flist = []

    for i, inf in enumerate(infs):
        with open(inf, 'r') as f:
            flines = f.readlines()
        data_2d = extract_numeric_block(flines, iline[0], iline[1]+1)
        ### dimension of input potential
        if inp_type == 2:
            flattend = np.array(data_2d).flatten()
        ngrid = len(list(flattend))
        v_g_flist.append(list(flattend))

    v_g_lists = np.array(v_g_flist)
    n_half = ngrid

    v_r2D = []
    print(f"size:v_g_lists {np.array(v_g_lists).shape}  length {ngrid}")
    for i in range(v_g_lists.shape[0]):
        #v_r = np.fft.ifft(np.fft.ifftshift(v_g_lists[i]))
        v_r = np.fft.ifft(v_g_lists[i])
        v_r2D.append(v_r)
    n_real = v_r.size

    k = np.linspace( 0, gmax, ngrid, endpoint=False)
    delta_k = gmax / (n_half - 1)
    delta_r = 2.0 * np.pi /(n_real * delta_k)
    r = np.arange(n_real) * delta_r

    y_real = np.array(v_r2D).T
    y_real = np.real(y_real)
    s = ""

    plot_list=[]

    for i, ncol in enumerate(y_real):
    #for i, ncol in enumerate(zip(v_g_flist, y_real)):
        print(f"{r[i]:>10.5f}", end=" ")
        s = f"{r[i]:>10.5f} "
        col_old = 0
        for col in ncol:
            print(f"{col:>15.7e}", end=" ")
            s += f"{col:12.7e} "
            if col_old:
                diff = col - col_old
            col_old = col
        if "diff" in locals(): 
            print(f"{diff:>15.7e}")
            s += f"{diff:>15.7e}\n"
        else:
            print()
            s += '\n'
        outf.write(s)
    #print(f"length: {len(flattend)}")
    return 0 
    

def main():

    parser = argparse.ArgumentParser(description="Extract pseudopotential from POTCAR ")
    parser.add_argument('-i','--inf',  nargs='*', help="several POTCARs")
    parser.add_argument('-o','--outf',  default='pseudo.dat', help="pseudopotential outfile")
    parser.add_argument('-t','--in_type', default=2, type=int, help="convert 2D array to 1D")
    parser.add_argument('-l','--iline', nargs=2,  type=int, default=[57, 256], help="pseudopot part in POTCAR")
    parser.add_argument('-l','--iline', nargs=2,  type=int, default=[57, 256], help="pseudopot part in POTCAR")
    parser.add_argument('-gm','--gmax', type=float, default=170.075338972103, help="to insert K-values in PSP at K-space")
    parser.add_argument('-m','--math', default='diff', help="calculate difference")
    #parser.add_argument('-r','--rspace', action='store_true', help="convert k-space to R-space")
    parser.add_argument('-j','--job',  choices=['xscale','fft', 'ifft'], help="job will make 1D file w. fft, ifft, etc")
    parser.add_argument('-u','--usage', action='store_true', help="show usage ")
    args = parser.parse_args()

    if args.usage:
        print("\
            \n\t\tpot_get_pseudo.py -i POTCAR2.H/0/POTCAR POTCAR2.H/0.1/POTCAR -l 57 256 > diff1.txt\
            \n\t\t(plot)\
            \n\t\t    gnu_1f_nc.sh diff1.txt 2 3 4\
            \n\t\t    gnu_1f_nc.sh diff1.txt 4\
            \n\t\t\
            \n")
        sys.exit(0)

    get_pseudopot(args.inf, args.in_type, args.iline, args.gmax, args.job, args.math, args.outf)

if __name__ == "__main__":
    main()
