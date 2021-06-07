#!/home/joonho/anaconda3/bin/python

import argparse
import sys
import re
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

def rearrange_coord_lines(lines, atoms, katoms):
    #print(katoms)
    ### from Null connect each line string to each index of atom kind
    newlines=[""]*len(katoms)
    for line, atom in zip(lines, atoms):
        ind = katoms.index(atom)
        newlines[ind]+=line
    #print(newlines[0])
    #print(newlines[1])
    ### connect each string (index) of newlines
    s=""
    for line in newlines:
        s+=line
    return s

def expand_atoms(atom_list, natom):
    atoms=[]
    times = natom / len(atom_list)
    atoms = atom_list * int(times)
    return atoms
    
### atom_list, natom || atom_file [ftype]
def poscar_rearrange(pos, atom_list, natom, atom_file, ftype):
    ### atom list is given repetively
    if atom_list:
        if natom == len(atom_list):
            atomlist=atom_list
        ### make full atomlist
        else:
            atomlist=expand_atoms(atom_list, natom)
    ### get atom list from bgf file 
    elif atom_file:
        if ftype == None:
            ftype = common.f_ext(atom_file)
        atomlist = get_atomlist4file(atom_file, ftype)
    else:
        print("one of -al and -af should be given")
        sys.exit(1)
    
    print(atomlist)
    ### rearrange poscar
    outf = open("POSCARnew", 'w')
    with open(pos, "r") as f:
        lines = f.readlines()
        i=0
        for line in lines:
            if i == 0:
                s = pos + "\n"
                outf.write(s)               # write() doesnot make "\n"
            elif i < 5:
                outf.write(line)
            ### arrange number of atoms with decreasing Atomic number
            elif i == 5:
                num_atom_kinds, katom_list = contract_atom_numbers(line, atomlist)
                s = "    " + "    ".join(katom_list) + "\n"
                outf.write(s)
                s = "    " + "    ".join(str(n_katom) for n_katom in num_atom_kinds) + "\n"
                outf.write(s)
            elif not re.search('\d', line):
                outf.write(line)
            else:
                break
            i+=1

    i_coord = i
    ### rearrange_coord_lines
    s = rearrange_coord_lines(lines[i_coord:], atomlist, katom_list)
    outf.write(s)

    return 0

def main():
    parser = argparse.ArgumentParser(description="rearange poscar atoms: (vmd,qchem)-poscar need to be rearranged")
    parser.add_argument('poscar', help="poscar to be modified")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-af','--atom_file', help="atom list in the order of original coordinate file ")
    group.add_argument('-al', '--atom_list', nargs='+', help="atom list is given directly and repetitively")
    parser.add_argument('-na', '--natom', type=int, help="if number of atoms are known")
    parser.add_argument('-ft','--file_type', help="original coordinate file type used with --atom_file")
    args = parser.parse_args()

    poscar_rearrange(args.poscar, args.atom_list, args.natom, args.atom_file, args.file_type)

if __name__ == "__main__":
    main()
