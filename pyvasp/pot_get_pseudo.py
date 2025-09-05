#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
import sys
import numpy as np
from common import dir_files
from my_data import extract_numeric_block
from libpotcar import parse_potcar 
from parsing import anylower, print_list2col
from mplot2D import mplot_nvector
#def extract_numeric_block(flines, i, j):
#    return [list(map(float, line.strip().split())) for line in flines[i:j]]
Lprint = 0

def get_pseudopot(infs, keys,  job, math, outfile, Lwrite_option):
    ### Make 2D lines

    if outfile:
        fwrite = []
        outf = open(outfile, 'w')


    pot_flist = []
    x = []
    labels=[]
    for ke in keys:
        for inf in infs:
            pot, buffer = parse_potcar(inf, kw=ke)
            ### dimension of input potential
            print(f"pot shape {np.array(pot).shape}")
            if isinstance(pot[0], list):
                x = [float(item) for item in np.array(pot)[:,0]]
                pot = [float(item) for item in np.array(pot)[:,1]]
            else:
                if ke == 'local':
                    gmax =  float(buffer)
                    if not x:
                        for i in range(len(pot)):
                            x.append(gmax/len(pot)*i)

            ngrid = len(pot)
            pot_flist.append(pot)
            label = inf + ke
            labels.append(label)
    
    
    #print(f"{pot_flist}")
    pot_arrs = np.array(pot_flist)
    print(f"shape pot_arrs: {pot_arrs.shape}")

    if job and 'fft' in job:
        v_r2D = []
    
        #print(f"size:v_g_lists {np.array(v_g_lists).shape}  length {ngrid}")

        for i in range(pot_arrs.shape[0]):
            v_g_sym = np.concatenate([pot_arrs[i], pot_arrs[i,-2:0:-1]])
            #v_r = np.fft.ifft(np.fft.ifftshift(v_g_lists[i]))
            v_r = np.fft.ifft(v_g_sym).real  # get real result
            N = len(v_r)
            v_r2D.append(v_r)

        #n_real = v_r.size

        #k = np.linspace( 0, gmax, ngrid, endpoint=False)
        delta_G = gmax / (N - 1)
        a = 2.0 * np.pi /(N * delta_G)
        r = np.linspace(0, a, N, endpoint=False)

        y_real = np.array(v_r2D).T
        y_real = np.real(y_real)

        #print(f"before plot r {r}")
        print(f"dimensions:: r {r.shape} v_r2D {np.array(v_r2D).shape}")
        mplot_nvector(r.tolist(), v_r2D, legend=labels)
        if outfile:
            for row in v_r2D.T.tolist(): #, y_real[0] ):
                s = f"{row[0]:>15.10f} {row[1]:>15.10f} \n"
                outf.write(s)

    else:
        print(f"before plot x dim {len(x)}")
        mplot_nvector(x, pot_arrs, legend=labels)
        if outfile:
            for xi, v_pots in zip(x, pot_arrs.T.tolist()): #, y_real[0] ):
                s = ""
                if Lwrite_option: s += f"{xi:>15.10f} "
                for ele in v_pots:
                    s += f"{ele:>15.10f} "
                s += "\n"
                outf.write(s)
    outf.close()
    #sys.exit(0)
    ### 0 for x, 1 for original, 2 for fft
    #if 1 in keys:
    #np_vglist = np.array(v_g_lists).T
    #plot_list=[]
    #plot_list.append(np_vglist[:,0])
    #print(f"plot_list {np.array(plot_list).shape}  length {ngrid}")

    ### plot list
    #if len(pot_arrs) == 1:
    #    ### whether x or not
    #print_list2col(pot_arrs[0], x=x, fname=outfile, Lfwrite=True)

    return 0 
    

def main():

    parser = argparse.ArgumentParser(description="Extract pseudopotential from POTCAR ")
    parser.add_argument('-i','--inf',  nargs='*', help="several POTCARs")
    parser.add_argument('-o','--outf',  default='pseudo.dat', help="pseudopotential outfile")
    parser.add_argument('-k','--key_list', nargs='*', default=['potcar'], choices=['potcar','siesta','alpha'],help="potfile type")
    parser.add_argument('-m','--math', default='diff', help="calculate difference")
    #parser.add_argument('-r','--rspace', action='store_true', help="convert k-space to R-space")
    parser.add_argument('-j','--job',  choices=['xscale','fft', 'ifft'], help="job will make 1D file w. fft, ifft, etc")
    parser.add_argument('-wo','--write_opt', action='store_true',  help="print x to outfile")
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

    get_pseudopot(args.inf, args.key_list, args.job, args.math, args.outf, args.write_opt)

if __name__ == "__main__":
    main()
