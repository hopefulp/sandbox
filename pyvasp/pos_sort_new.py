#!/home/joonho/anaconda3/bin/python

import argparse
import sys
import re
import os
import numpy as np

import common
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

def rearrange_coord_lines(lines, atom_names, natoms):
    #print(katoms)
    ### from Null connect each line string to each index of atom kind
    coord_dict={}
    for aname, natom in zip(atom_names, natoms):
        if aname not in coord_dict.keys():
            coord_dict[aname]=[]
        coord_dict[aname].extend(lines[:natom])
        del lines[:natom]
    s=''
    for atom in coord_dict.keys():
        for line in coord_dict[atom]:
            print(line, end='')
            s += line
    return s

def expand_atoms(atom_list, natom):
    ''' make it fold for atoms '''
    atoms=[]
    times = natom / len(atom_list)
    atoms = atom_list * int(times)
    return atoms

def cal_num_atoms(alist, nalist):
    mydic = {}
    print(alist, nalist)
    for atom, natom in zip(alist, nalist):
        if atom in mydic.keys():
            mydic[atom] = mydic[atom] + natom # += is not working
        else:
            mydic[atom] = natom
    return mydic            

### Just read poscar and rearrange: atom_list, natom || atom_file [ftype]
def poscar_rearrange(pos, atom_list, natom, atom_file, ftype, suff, rename):
    ### get atom list from bgf file 
    if atom_file:
        if ftype == None:
            ftype = common.f_ext(atom_file)
        atomlist = get_atomlist4file(atom_file, ftype)
    #else:
    #    print("one of -al and -af should be given")
    #    sys.exit(1)
    
    #print(atomlist)
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
                if not atom_list: 
                    if len(line.split()) >= 2:
                        alist_all = line.strip().split()
                        atoms = []
                        for atom in alist_all:
                            if atom not in atoms:
                                atoms.append(atom)
                    el_st = ''.join(atoms)
                    s = el_st + "   from " + pos + "\n"
                    #print(s)
                else:
                    s = 'atomlist from outside'
                outf.write(s)               # write() doesnot make "\n"
            ### copy cell parameters
            elif i < 5:
                outf.write(line)
            ### filter unique atom species in the order in POSCAR
            elif i == 5 and not any(s.isdigit() for s in line) :
                ### keep atom list with index
                atoms_all = line.strip().split()
                s = "  " + "  ".join(atom_list) + "\n"
                outf.write(s)
            ### calculate sum of atom-kinds
            elif i == 6 and any(d.isdigit() for d in line):
                natom_all = list(map(int, line.strip().split()))
                number_dic = cal_num_atoms(atom_all, natom_all)
                print(number_dic)
                s=''
                number_list=[]
                for atom in atoms:
                    s += f"  {number_dic[atom]} "
                    number_list.append(number_dic[atom])
                ntatom = sum(number_list)
                s += "\n"
                outf.write(s)
            elif not re.search('\d', line):
                outf.write(line)
            ### sort coordinates
            else:
                coordinates.append(line)
    print(f"{ntatom} {len(coordinates)}")
    ### rearrange_coord_lines
    s = rearrange_coord_lines(coordinates, atom_all, natom_all)
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
    parser.add_argument('-al', '--atom_list', nargs='+', help="atom list is given directly and repetitively")
    parser.add_argument('-na', '--natom', type=int, help="if number of atoms are known")
    parser.add_argument('-ft','--file_type', help="original coordinate file type used with --atom_file")
    parser.add_argument('-suf', '--suffix', default='sort', help='filename suffix')
    parser.add_argument('-mv', '--rename', action='store_true', help='change in the original filename')
    args = parser.parse_args()

    poscar_rearrange(args.poscar, args.atom_list, args.natom, args.atom_file, args.file_type, args.suffix, args.rename)

if __name__ == "__main__":
    main()
