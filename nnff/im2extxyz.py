#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os
import numpy as np
from common import *
from _math import minmax
import molecules
from my_chem import C100kj2ev

### scratch keywords for AMP from IM/MM input file of the form of .fin
kw_ene = "energy"
kw_coord = "cartesian"

def extract_data(fin, atom_list, ny_bar, dmining):

    #print(fin)
    #print(os.path.splitext(fin))
    #print(molecules.lattice_properties)
    natoms = len(molecules.mol2atoms[atom_list])
    #print(natoms)              # check natoms from molecules.py
    if ny_bar == 1:
        comment_pre = molecules.lattice_properties_e + "energy="
    elif ny_bar == 2:
        comment_pre = molecules.lattice_properties_ef + "energy="
    f_extxyz = os.path.splitext(fin)[0] + ".extxyz"
    fout=open(f_extxyz, 'w')

    #print("natoms = %d" % len(molecules.atom_list))

    ### Machine Learning params
    y_max=0.0
    y_min=0.0

    ### loop control params
    tag="find_ene"
    iconfig=0
    with open(fin, 'r') as f:
        lines=f.readlines()
        i=0
        for line in lines:                          # line has "\n"
            # find energy tag then find cartesian coordinate alternatively
            if tag == "find_ene":
                if re.search(kw_ene, line):
                    tag="find_coord"
                    print(natoms, end="\n")
                    natoms_line = str(natoms) + "\n"
                    fout.write(natoms_line)
                    # complete comment line
                    y_value =  float(lines[i+1]) * C100kj2ev
                    if dmining:
                        y_value = y_value * float(10**dmining)
                    comment_line = comment_pre + "%-10.5f"%y_value +"\n"   # add energy to line
                    if y_min==0.0 and y_max==0.0:
                        y_min=y_max=y_value
                    else:
                        y_min, y_max = minmax(y_value, y_min, y_max)
                    print(comment_line,end="")
                    fout.write(comment_line)
                    iconfig+=1
                    i+=1
                    continue
                else:
                    i+=1
                    continue
            # find cartesian tag        
            elif tag == "find_coord":
                if re.search(kw_coord, line):
                    # write natom lines: treat block in advance
                    j=1
                    for ia in molecules.mol2atoms[atom_list]:
                        Z = molecules.atomic_number[ia]
                        atom_line = "%-5s"%(ia)+lines[i+j].strip()+"%5d"%(Z)+" %5.1f %5.1f %5.1f"%(0.0,0.0,0.0) + "\n"
                        print(atom_line,end="")
                        fout.write(atom_line)
                        j+=1
                    tag = "find_ene"
                    i+=1
                    continue
                else:
                    i+=1
                    continue

    fout.close()
    print("sampled {} configurations".format(iconfig))
    print("min = {}, max={}".format(y_min, y_max))

    return

def main():
    parser = argparse.ArgumentParser(description='make extxyz from IM input file (.fin)')
    parser.add_argument('fin', help='IM/MM input file of .fin')
    parser.add_argument('-a', '--atoms', default='ethylene', choices=['ethylene','Diss_CHO'], help='atom list of the molecule')
    parser.add_argument('-y', '--ny_bar', default=1, type=int, choices=[1,2,3], help='machine learning target={1:ene,2:ene & force,3:ene & force & hessian}')
    parser.add_argument('-d', '--dmining', type=int, choices=[1,2,3,4], help='data mining: value*10^n')
    args = parser.parse_args()
    extract_data(args.fin, args.atoms,args.ny_bar,args.dmining) 

if __name__ == '__main__':
    main()
