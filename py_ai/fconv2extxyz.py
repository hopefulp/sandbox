#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os
import numpy as np
from common import *
from _math import minmax
import chem_space
from my_chem import C100kj2ev, hr2ev

### scratch keywords for AMP from IM/MM input file of the form of .fin
kw_ene = "energy"
kw_coord = "cartesian"

def from_xyz(fin, atom_list, y_bar_type, dmining):
    chem_sys = chem_space.H2O()

    ext_str = chem_sys.extxyz_str+" energy="
    #print (chem_sys.extxyz_str)
    f_extxyz = fin.split('.')[0]+".extxyz"
    fout=open(f_extxyz, 'w')

    ### Machine Learning params
    y_max=0.0
    y_min=0.0

    with open(fin, 'r') as f:
        lines=f.readlines()
        i=0
        for line in lines:
            if i==0 :
                natom = int(line.strip())
                fout.write(line)
                #print(line, end='')
                i+=1
                continue
            if i%(natom+2)==0:
                fout.write(line)
                #print(line, end='')
                i+=1
                continue
            if re.search("E", line):
                y_value = float(line.strip().split()[-1]) * hr2ev
                pline = ext_str + str(y_value) + "\n"
                if y_min==0.0 and y_max==0.0:
                    y_min=y_max=y_value
                else:
                    y_min, y_max = minmax(y_value, y_min, y_max)
                fout.write(pline)
                #print(pline,end='')
                i+=1
                continue
            if 2 <= i%(natom+2):
                ele = line.strip().split()
                atom_line = "{:2} {:10.6} {:10.6} {:10.6}".format(ele[0],float(ele[1]),float(ele[2]),float(ele[3]))
                Z = chem_space.atomic_number[ele[0]]
                atom_line += "%5d"%(Z)+" %5.1f %5.1f %5.1f"%(0.0,0.0,0.0) + "\n"
                #print(atom_line, end='')
                fout.write(atom_line)
                i+=1
                continue

def from_fin(fin, atom_list, y_bar_type, dmining):
    natoms = len(chem_space.mol2atoms[atom_list])
    if y_bar_type == 1:
        comment_pre = chem_space.lattice_properties_e + "energy="
    elif y_bar_type == 2:
        comment_pre = chem_space.lattice_properties_ef + "energy="
    f_extxyz = os.path.splitext(fin)[0] + ".extxyz"
    fout=open(f_extxyz, 'w')
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
                    for ia in chem_space.mol2atoms[atom_list]:
                        Z = chem_space.atomic_number[ia]
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
    parser = argparse.ArgumentParser(description='obtain extxyz from IM input file (.fin)')
    parser.add_argument('fin', help='IM/MM input file of .fin')
    parser.add_argument('-a', '--atoms', default='ethylene', choices=['ethylene','Diss_CHO'], help='atom list of the molecule')
    parser.add_argument('-y', '--ny_bar_type', default=2, type=int, choices=[1,2,3], help='machine learning target={1:ene,2:ene & force,3:ene & force & hessian}')
    parser.add_argument('-d', '--dmining', type=int, choices=[1,2,3,4], help='data mining: value*10^n')
    args = parser.parse_args()

    f_ext = args.fin.split('.')[-1]
    if f_ext == 'fin':
        from_fin(args.fin, args.atoms, args.ny_bar_type, args.dmining)
    elif f_ext == 'xyz':
        from_xyz(args.fin, args.atoms, args.ny_bar_type, args.dmining)
if __name__ == '__main__':
    main()
