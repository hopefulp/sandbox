#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os
from common import fname_decom

def formatt(el1, el2):
    st = f'{el1:>10.3f}{el2:12.4f}\n'
    return st

def math_file(fname, math, math_value, xcol, ycol, xopt, xpivot):
    '''
    xopt, xpivot: [str, float] or [float, float] or [int, int]
    xopt    gt, lt, float, int
    xpivot  0.5,    float, int
    '''
    #print(f'xcol {xcol}; ycol {ycol}; xopt {xopt}; xpivot {xpivot}')
    fn, fsuff = fname_decom(fname)
    ofname = f'{fn}_new{int(math_value)}.{fsuff}'
    fout=open(ofname, 'w')
    new_lines=[]
    with open(fname, 'r') as f:
        ### analyze file
        lines = f.readlines()
        nline = len(lines)
        ### get range
        for i, line in enumerate(lines):
            eles = list(map(float, line.strip().split()))
            if xopt == 'gt':
                if eles[xcol] < xpivot:
                    newline = line      # f.write() has no newline
                else:
                    if math == 'prod':
                        value = eles[ycol] * math_value
                    elif math == 'div':
                        value = eles[ycol] / math_value
                    newline = formatt(eles[xcol], value)
            new_lines.append(newline)
            ### directly write to new file
            fout.write(newline)
    fout.close()            
    return 0

def main():
    parser = argparse.ArgumentParser(description='To get rid of abnormal DOS at start energy')
    parser.add_argument('file',  help='input file')
    parser.add_argument('-m', '--math', default='prod', choices=['prod','add','sub','div'], help='math operator on column values')
    parser.add_argument('-v', '--math_value', type=float, default=5.0, help='value for operator')
    parser.add_argument('-x', '--xcol', default=0, type=int,  help='column index')
    parser.add_argument('-y', '--ycol', default=1, type=int,  help='column index')
    gregion = parser.add_mutually_exclusive_group()
    gregion.add_argument('-hlt', '--half_less_than', help='half-less than in 0-th column indices')
    gregion.add_argument('-hgt', '--half_greater_than', help='half-greater than in 0-th column indices')
    gregion.add_argument('-xgt', '--x_greater_than', type=float, help='greater than x in 0th-column indices')
    gregion.add_argument('-xlt', '--x_less_than', type=float, help='less than x in 0th-column indices')
    gregion.add_argument('-xij', '--xij_index', nargs=2, type=float, help='x region')
    args = parser.parse_args()

    if args.half_less_than:
        xopt = 'lt'
        xpivot = 0.5
    elif args.half_greater_than:
        xopt = 'gt'
        xpivot = 0.5
    elif args.x_less_than:
        xopt = 'lt'
        xpivot = args.x_less_than
    elif args.x_greater_than:
        xopt = 'gt'
        xpivot = args.x_greater_than
    elif args.xij_index:
        xopt = args.xij_index[0]
        xpivot = args.xij_index[1]
    
    math_file(args.file, args.math, args.math_value, args.xcol, args.ycol, xopt, xpivot) 

if __name__ == '__main__':
    main()
