#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os
from _extract_line import extract_values
from mod_poscar import parse_poscar

def find_igroup(k, nlist):
    accum = 0
    for i, n in enumerate(nlist):
        accum += n
        if k <= accum:
            return i


def cal_pcharge(fname, atoms, natoms, zvals):

    if not zvals:
        zvals = extract_values('POTCAR', 'ZVAL')
    print(f"zval = {zvals}")
    tatom_list, tnatom_list = parse_poscar('POSCAR', 'alist')
    if not atoms:
        atoms = tatom_list
        natoms = tnatom_list
    print(f"atoms {atoms} natoms {natoms}")
    ofname = "pcharge_bader.dat"
    
    fout=open(ofname, 'w')
    with open(fname, 'r') as f:
        lines = f.readlines()

        for i, line in enumerate(lines):          # line has "\n"
            #print line,         # print writes its own "\n"
            eles = line.strip().split()
            if eles[0].isdigit():
                ind = int(eles[0])
                igroup = find_igroup(ind, tnatom_list)
                pchg = float(zvals[igroup]) - float(eles[4])
                s = f" {ind:>2} {float(zvals[igroup]):5.2f} {pchg:6.3f}"
                print(s)
                fout.write(s+"\n")    # write does not write "\n"
    return 0

def main():
    parser = argparse.ArgumentParser(description='Partial charge w.r.t. ZVAL in POTCAR')
    parser.add_argument('file', nargs='?', default='ACF.dat', help='read ACF.dat')
    parser.add_argument('-a', '--atom_list', nargs='*', help='atom list')
    parser.add_argument('-n', '--natom_list', nargs='*', help='natom list')
    parser.add_argument('-z', '--zval_list', nargs='*', help='zval list')

    args = parser.parse_args()
    
    cal_pcharge(args.file, args.atom_list, args.natom_list, args.zval_list) 

if __name__ == '__main__':
    main()
