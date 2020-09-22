#!/gpfs/home/joonho/anaconda3/bin/python

import argparse
import os
import re
import numpy as np
#import itertools
import numpy as np
import my_stat
from common import dir_files

file_ene = 'test_energy.txt'
file_force = 'test_force_img_000.dat'

def get_inf_name(ef):
    if not ef:
        print("if not input file w. -inf, input -e|-f|-ef ")
        sys.exit(1)
    if re.search('e', ef):
        if os.path.isfile(file_ene):
            f = file_ene
    if re.search('f', ef):
        if os.path.isfile(file_force):
            f = file_force
    return f            

def run_stat(inf, stat_ef, scale):
    '''
    scale in error
        energy: scale = power(10,6) as np.square(meV)
        force : scale = 1           as eV/A
    '''
    if not inf:
        inf = get_inf_name(stat_ef)

    data_line=[]
    x  = []
    y  = []
    if 'e' in stat_ef:
        with open(inf, 'r') as f:
            lines = f.readlines()
            #for line in itertools.islice(f, 1, None):
            for line in lines:
                if re.search('[a-zA-Z]', line):
                    continue
                ele = line.strip().split()
                x.append(float(ele[0]))
                y.append(float(ele[1]))
        p = my_stat.stat_2col(x, y, scale=np.power(10,6))
        ffmt = '9.4f'
        sfmt  = '^8'
        print(f"{'MSE':{sfmt}}{p.mse:{ffmt}}\n{'BIAS_SQ':{sfmt}}{p.bias_sq:{ffmt}}\
                \n{'Var(X-Y)':{sfmt}}{p.varx_y:{ffmt}}    {'Var(true)':{sfmt}}{p.varx:{ffmt}}    {'Var(hypo)':{sfmt}}{p.vary:{ffmt}}    {'Cov(X,Y)':{sfmt}}{p.cov:{ffmt}}    {'r_corr':{sfmt}}{p.r_corr:10.4f}")

    elif 'f' in stat_ef:
        with open(inf, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if re.search('[a-zA-Z]', line):
                    continue
                ele = line.strip().split()
                ### combine f_x, f_y, f_z in x[]
                for i in range(3):
                    x.append(float(ele[2*i]))             # force_x x.append(float(ele[0]))
                    y.append(float(ele[2*i+1]))
        p = my_stat.stat_2col(x, y, scale=1)
        ffmt = '>12.4f'
        sfmt  = '^8'
        print(f"{'MSE':{sfmt}}{p.mse:{ffmt}}\n{'BIAS_SQ':{sfmt}}{p.bias_sq:{ffmt}}\
                \n{'Var(X-Y)':{sfmt}}{p.varx_y:{ffmt}}    {'Var(true)':{sfmt}}{p.varx:{ffmt}}    {'Var(hypo)':{sfmt}}{p.vary:{ffmt}}    {'Cov(X,Y)':{sfmt}}{p.cov:{ffmt}}    {'r_corr':{sfmt}}{p.r_corr:10.4f}")

    return 0 
    

def main():

    parser = argparse.ArgumentParser(description="analyse data: error, variance  ")
    parser.add_argument('-inf', '--inputfile',  help=" input file ")
    parser.add_argument('-s', '--scale', help="use E-6 for energy(meV)| 1 for force(eV/A)")
    group_cal = parser.add_mutually_exclusive_group()
    group_cal.add_argument('-e', '--energy', action='store_true', help='calculate statistics for energy')
    group_cal.add_argument('-f', '--force',  action='store_true', help='calculate statistics for force')
    #group_cal.add_argument('-ef', '--ene_force', action='store_true', help='calculate statistics for energy and force')
    args = parser.parse_args()

    st = ''
    if args.energy:
        st = 'e'
    elif args.force:
        st = 'f'

    run_stat(args.inputfile, st, args.scale)

if __name__ == "__main__":
    main()
