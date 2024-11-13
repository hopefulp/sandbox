#!/home/joonho/anaconda3/bin/python
import argparse
import re
import os
import sys
from common import list2dict, whereami
from libincar import modify_incar_bykv

def change_incar_bydic(incar, kws, outf='incar'):
    '''keeping file, change line'''
    modify_incar_bydic(incar, kws, outf='incar')
    return 0

def change_incar_byjob(incar, job, outf=None):
    if job == 'cont':
        kw = {'ISTART': 1, 'ICHARG': 0}
    else:
        print("Input job is required")
    
    ### if outfile == INCAR, change old file
    if outf == 'INCAR':
        os.system("cp INCAR INCAR_o")
    #change_incar_bydic(incar, kw, outf=outf)
    modify_incar_bydic(incar, kw, outf=outf)
    return 0


def main():
    parser = argparse.ArgumentParser(description='change INCAR parameter')
    parser.add_argument('inf', help='input incar file or directory')
    parser.add_argument('-j', '--job', choices=["dos","band","pchg","chg","md","cont","ini","zpe","mol","wav",'vdw','noD','opt','copt','mag','kisti'], help='job for VASP')
    parser.add_argument('-kv', '--keyvalue', help='input kw to change') 
    parser.add_argument('-o', '--output',  help='change the output filename')
    args = parser.parse_args()

    if os.path.isfile(args.inf):
        f = args.inf
    elif os.path.isdir(args.inf):
        f = args.inf + "/INCAR"

    if args.output:
        outf = args.output
    else:
        outf = 'INCAR_n'
        
    if args.job:
        change_incar_byjob(f, args.job, outf=outf)
    elif args.keyvalue:
        ### change list2dict
        dic = list2dict(args.keyvalue)
        change_incar_bydic(f, dic, outf=outf)
    return 0

if __name__ == "__main__":
    main()

