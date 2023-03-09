#!/home/joonho/anaconda3/bin/python

import argparse
import sys
import re
import os
import numpy as np

from common import whereami
import chem_space as cs

def get_atomlist4bgf(fname):
    keyword1 = "HETATM"
    keyword2 = "ATOM"
    atomlist=[]
    #print(fname)
    with open(fname, 'r') as f:
        lines = f.readlines()
        flag='NO'
        
        for line in lines:
            if not (re.match(keyword1, line) or re.match(keyword2,line)):
                continue
            else:
                l_ele = line.split()
                index_d = re.search('\d',l_ele[2])
                #print(f"{l_ele[1]} {l_ele[2]} {index_d.start()}")
                atomlist.append(l_ele[2][:index_d.start()])
    
    print(("natom = ", len(atomlist)))
    print(atomlist)
    return atomlist

def get_atomlist4file(fname, ftype):
    if ftype == "bgf":
        alist = get_atomlist4bgf(fname)
    else:
        print(f"{ftype} file type is not ready")
        sys.exit(1)
    return alist

def get_atom_kinds(atoms):
    atom_k = []
    for atom in atoms:
        if not atom in atom_k:
            #print(atom)
            atom_k.append(atom)
    atom_k = cs.arrange_atom_list(atom_k)            
    print(atom_k)
    return atom_k

### return number of atoms as for atom kinds and atom kinds itself
def contract_atom_numbers(numbers, atomlist):
    atom_kinds = get_atom_kinds(atomlist)
    num_atom_kinds=[0]*len(atom_kinds)
    for atom in atomlist:
        ind = atom_kinds.index(atom)
        num_atom_kinds[ind]+=1
    print(num_atom_kinds)
    return num_atom_kinds, atom_kinds

def sortz_atom_coord(lines, atom):
    ''' for a same atom group '''
    line2d=[]
    for line in lines:
        coords = line.strip().split()   # sometimes including T T F
        line2d.append(coords)
    #print(f"{line2d} in {whereami()}")
    new_coord2d = sorted(line2d, key=lambda l:float(l[2]))
    ### transform into 1d
    s=[]
    for line_ele in new_coord2d:
        ### use join to make a string regardless of number of data
        ### sorted POSCAR format
        space = "     "
        icoord = list(map(float, line_ele))
        st = f"  {icoord[0]:10.5f} {icoord[1]:10.5f} {icoord[2]:10.5f}   {atom}\n"
        s.append(st)
    return s

def rearrange_coord_lines(lines, atom_klist, atom_names_orig, natom_orig, sort_z):
    #print(katoms)
    ### from Null connect each line string to each index of atom kind
    coord_dict={}
    ### make dict[atom_name]=[coordinate_string, ... ]
    for aname, natom in zip(atom_names_orig, natom_orig):
        if aname not in coord_dict.keys():
            coord_dict[aname]=[]
        coord_dict[aname].extend(lines[:natom])
        del lines[:natom]
    s=''
    for atom in atom_klist:
        if sort_z:
            print(f"z sort in {whereami()}")
            lines = sortz_atom_coord(coord_dict[atom], atom=atom)
        else:
            lines = coord_dict[atom]
        for line in lines:
            print(line, end='')
            s += line
    return s

def expand_atoms(atom_list, natom):
    ''' make it fold for atoms '''
    atoms=[]
    times = natom / len(atom_list)
    atoms = atom_list * int(times)
    return atoms

### Below: used in poscar_sort()

def cal_num_atoms_orig(alist, nalist):
    ''' return
            dict of dict[atoms_species]=natom for original poscar
    '''
    mydic = {}
    print(alist, nalist)
    for atom, natom in zip(alist, nalist):
        if atom in mydic.keys():
            mydic[atom] = mydic[atom] + natom # += is not working
        else:
            mydic[atom] = natom
    return mydic            

### calculate sum of atom-kinds
def obtain_sort_natom_list(line, katom_orig, atom_klist):
    natom_orig = list(map(int, line.strip().split()))
    number_dic = cal_num_atoms_orig(katom_orig, natom_orig)
    print(number_dic)
    s=''
    sort_natom_list=[]
    for atom in atom_klist:
        s += f"  {number_dic[atom]} "
        sort_natom_list.append(number_dic[atom])
    ntatom = sum(sort_natom_list)
    s += "\n"
    return s, ntatom, natom_orig


def get_atomlist(line):
    alist_all = line.strip().split()
    atoms = []
    for atom in alist_all: # also needs to check whether it is atom symbol
        if atom not in atoms:
            atoms.append(atom)
    return atoms, alist_all  

### Just read poscar and rearrange: atom_list, natom || atom_file [ftype]
### POSCAR sort for the same atoms
def poscar_sort(pos, atom_klist, natom, atom_file, ftype, suff, rename, sort_z):
    ### atom list is given repetively --> atom kinds list
    ### rearrange poscar
    ofile = pos+suff
    outf = open(ofile, 'w')
    ### extract coordinates and sort
    coordinates=[]
    with open(pos, "r") as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if i == 0:
                ### if the 0-th line is atom list, get atom list
                if len(line.split()) >= 2:
                    atom_sort, katom_orig = get_atomlist(line)
                if not atom_klist:
                    atom_klist = atom_sort
                el_st = ''.join(atom_klist)
                s = el_st + "   from " + pos + "\n"
                print(s)
                outf.write(s)               # write() doesnot make "\n"
            ### copy cell parameters
            elif i < 5:
                outf.write(line)

            ### filter unique atom species in the order in POSCAR
            ### there might be atom species or not
            elif i == 5:
                ### if atom symbols line exists
                if not any(s.isdigit() for s in line):
                    ### keep atom list with index
                    atom_sort, katom_orig = get_atomlist(line)
                    if not atom_klist:
                        atom_klist = atom_sort
                    s = "  " + "  ".join(atom_klist) + "\n"
                    print(s)
                    outf.write(s)
                    tag_atomsort = False ##--> read natoms line
                ### if not, calculate natoms list
                else:
                    s, ntatom, natom_orig = obtain_sort_natom_list(line, katom_orig, atom_klist)
                    outf.write(s)
                    tag_atomsort = True
            ### in case: at i==6, natoms exist
            elif not tag_atomsort:
                s, ntatom, natom_orig = obtain_sort_natom_list(line, katom_orig, atom_klist)
                outf.write(s)
                tag_atomsort = True

            ### write cartesian or direct if not digit
            elif re.match('[SDC]', line, re.I):     # write and skip if Selective|Direct|Cartesion
                print(f"line: {line} in {whereami()}")
                outf.write(line)
            ### sort coordinates
            elif len(coordinates) < ntatom:
                coordinates.append(line)
            else:
                ### if POSCAR does not end with coordinates
                break
    print(f"total {ntatom} atoms with {len(coordinates)}-coordinates in original")
    ### rearrange_coord_lines
    s = rearrange_coord_lines(coordinates, atom_klist, katom_orig, natom_orig, sort_z)
    outf.write(s)
    if rename:
        os.system(f"mv {ofile} {pos}")
        print(f"{pos} was changed")
    else:
        print(f"write to {ofile}")

    return 0

def main():
    parser = argparse.ArgumentParser(description="rearange poscar atoms: (vmd,qchem)-poscar need to be rearranged")
    parser.add_argument('poscar', help="poscar to be modified")
    parser.add_argument('-af','--atom_file', help="atom list in the order of original coordinate file ")
    #parser.add_argument('-al', '--atom_list', nargs='+', help="atom list is given directly and repetitively")
    parser.add_argument('-al', '--atom_klist', nargs='+', help="atom kind list")
    parser.add_argument('-na', '--natom', type=int, help="if number of atoms are known")
    parser.add_argument('-ft','--file_type', help="original coordinate file type used with --atom_file")
    parser.add_argument('-suf', '--suffix', default='s', help='filename suffix')
    parser.add_argument('-mv', '--rename', action='store_true', help='change in the original filename')
    parser.add_argument('-z', '--sortz', action='store_true', help='increasing order in z-axis')
    args = parser.parse_args()

    poscar_sort(args.poscar, args.atom_klist, args.natom, args.atom_file, args.file_type, args.suffix, args.rename, args.sortz)

if __name__ == "__main__":
    main()
