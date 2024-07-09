#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
import sys

def make_body(line, atoms):
    body_list=[]
    ele = line.strip().split()
    body_list.append(f"{float(ele[0]):8.3f}\n")
    i=0
    for atom in atoms:
        body_list.append(f"{atom:<2s}{float(ele[i+1]):15.10f}{float(ele[i+2]):15.10f}{float(ele[i+3]):15.10f}\n")
        i+=3
    return body_list


def convert(dname, atoms):
    natoms = len(atoms)
    wd = os.getcwd()
    if dname:
        wd += '/' + dname 
    inf = wd + '/NucCarts'
    if not os.path.isfile(inf):
        print(f"No NucCarts in {wd}")
        sys.exit(1)
    outf = wd + '/View.xyz'
    of = open(outf, 'w')
    with open(inf, 'r') as f:
        lines = f.readlines()
        i=0
        for line in lines:
            if i==0:
                i+=1
                continue
            else:
                of.write(f"{natoms}\n")
                for aline in make_body(line, atoms):
                    of.write(aline) 
    return 0

def main():

    parser = argparse.ArgumentParser(description="convert NucCarts to xyz")
    parser.add_argument('-d', '--dirname', help="directory for file location ")
    parser.add_argument('-a','--atoms', nargs='*', help="atom series ")
    args = parser.parse_args()

    if args.atoms==None:
        print("input atom list")
        sys.exit(0)

    convert(args.dirname, args.atoms)

if __name__ == "__main__":
    main()

