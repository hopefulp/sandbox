#!/home/joonho/anaconda3/bin/python
'''
    made by J. Park 2022.01.03 runs for POSCAR, KPOINTS, PROCAR
'''
import argparse
import os
import re
import sys
import numpy as np
from common         import dir_files
from mod_poscar     import read_poscar
from mod_kpoints    import read_kpoints

def get_iatoms_in_group(zmin, zmax, coord, job):
    ind=[]
    for i, xyzs in enumerate(coord):
        line_ele = xyzs.strip().split()
        print(f"{line_ele}")
        if job == 'bulk':
            if zmin < float(line_ele[2]) and float(line_ele[2]) < zmax:
                ind.append(i)
        elif job == 'surf':
            if float(line_ele[2]) < zmin or zmax < float(line_ele[2]):
                ind.append(i)
    return ind

def obtain_atomlist(zminmax, poscar, atom_species, job):
    '''
    input
        read poscar
    return
        atom list inbetween zmin & zmax
        principal axes    
    '''
    if len(zminmax) == 2:
        zmin, zmax = (zminmax[0], zminmax[1])
    else:
        print("z-coord error: {zminmax}, input two z-values with -z ")
        return 1
    atom_list, natom_list, coord, direct_z = read_poscar('atomcoord', poscar)
    ### in case direct coordinates
    if direct_z != 1.0:
        zmin /= direct_z
        zmax /= direct_z
        print(f"zmin, max = {zmin} {zmax}")
    iatom=0
    ind_select=[]
    print(f"{atom_list} {natom_list} {len(coord)}")
    for i, atom in enumerate(atom_list):
        if atom in atom_species:
            ind_group = get_iatoms_in_group(zmin, zmax, coord[iatom:iatom+natom_list[i]], job)
            ind_select.extend([x+iatom for x in ind_group])
        iatom += natom_list[i]
            
    print(f"indices {ind_select} total {len(ind_select)}")
    return ind_select

def cal_kspacing(paxes2d):
    '''
    input
        principal axis to calculate spacing
        read KPOINTS
    output
        npoints for each kline
        number of kpoint lines
        pi/a following each k-line
    '''

    npoints, nkline = read_kpoints('KPOINTS')
    print(f'npoints and nkline: {npoints} {nkline}')
    kdir = ['y', 'x', 'y', 'x']
    pi_ax = np.pi/paxes2d[0][0]
    pi_ay = np.pi/paxes2d[1][1]
    
    pi_a = []
    for kd in kdir:
        if kd == 'x':
            pi_a.append(pi_ax)
        elif kd == 'y':
            pi_a.append(pi_ay)
    spacing = [ x/(npoints -1) for x in pi_a ]
    return npoints, nkline, spacing
        

def analyze_procar(atom_list, poscar, kpoints, procar):

    ### make kpoints distance
    p_axes = read_poscar('paxes', poscar)
    ##  For the given paths: S - X - G - Y - S: y-dir, x-dir, y-dir, x-dir
    npoints, nkline, spacing = cal_kspacing(p_axes)
    print(f"npoints {npoints}, nkline {nkline}, spacing {spacing}") # check cal_kspacing
    nheadline = 2 
    with open('PROCAR', 'r') as f:
        lines = f.readlines()       # load all the lines in PROCAR
        state = lines[1].split()
        nkpoints    = int(state[3])
        nbands  = int(state[7])
        nions   = int(state[11])
        i  = 0
        i += nheadline
        for k in range(nkpoints):
            ### save k
            i += 2
            for b in range(nbands):
                ### save b
                i += 2
                ### table label
                i += 1
                print(f"{lines[i]} for label in k {k} band {b}")
                for ia in range(nions):
                    paband = lines[i+ia].strip().split()
                    if paband[0] in atom_list:
                        print(f"{paband} in atoms loop")
                ### skip total line at each k, b and empty line
                i += nions + 2
            ### for k-block spacing
            i += 1 
    return 0 
    

def main():

    parser = argparse.ArgumentParser(description="display Usage for /mymplot  ")
    parser.add_argument('procar', nargs='?', default="PROCAR", help="input PROCAR")
    parser.add_argument('poscar', nargs='?', default="POSCAR", help="input POSCAR")
    parser.add_argument('kpoints', nargs='?', default="KPOINTS", help="input KPOINTS file for band")
    parser.add_argument('-j','--job', default='bulk',choices=['bulk','surf','ads'], help="selective band structure")
    atoms = parser.add_mutually_exclusive_group()
    atoms.add_argument('-al','--atom_list', nargs='*',  help="list atoms ")
    atoms.add_argument('-z', '--zintv', nargs=2, type=float, help='atoms between zmax * zmin')
    parser.add_argument('-as', '--atom_species', nargs='+', default=['O'], help='specify atom species')
    args = parser.parse_args()

    ### obtain atom list
    if args.atom_list:
        alist = args.atom_list
    else:
        if not args.zintv:
            print("Error: input z values for zmin and zmax")
            sys.exit(1)
        else:
            alist = obtain_atomlist(args.zintv, args.poscar, args.atom_species, args.job)
    analyze_procar(alist, args.poscar, args.kpoints, args.procar)

if __name__ == "__main__":
    main()
