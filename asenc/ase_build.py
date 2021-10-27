#!/home/joonho/anaconda3/bin/python

import argparse
from ase import Atoms, Atom, units
import ase.build as bd
from ase.visualize import view
from ase.io.formats import read, iread, write, string2index
from random import uniform
import numpy as np
import sys
import re
from ase.cluster import *
from common import list2str

class Usage:
    pass

metals = {'Pt': 'fcc'}

def make_random_coordinates(natoms, alattice):
    atoms_coord=[]
    for i in range(natoms):
        a_coord=[]
        for j in range(3):
            coord = uniform(0, alattice)
            a_coord.append(coord)
        atoms_coord.append(a_coord)
    return atoms_coord

def make_ordered_coordinates(natoms, alattice):
    atoms_coord=[]
    nx = int(np.pow(natoms, 1/3))

    #for i in range(nx):
            
def build_structure(dim, structure, name, size, vac, akind, status=None,  atom=10, latt=50):
    if len(size) == 1:
        size0 = size[0]
        size = (size0, size0, size0) 
    if len(vac) == 1:
        vac0 = vac[0]
        vac  = (vac0, vac0, vac0)
    if dim == 0:
        if re.match('icosa', structure):
            atoms = Icosahedron(akind, size0)
        elif re.match('deca', structure):
            atoms = Decahedron(akind, p=size[0], q=size[1], r=0)
            print(f"Decahedron with p {size[0]} q {size[1]} and r=0")
        elif re.match('octa', structure):
            atoms = Octahedron(akind, length=size0)
        else:
            print(f"Input cluster structure with -s {structure} such as: icosa,deca ")
            sys.exit(1)
        atoms.set_cell(vac)
        tr2d = [ [ vac0/2. for i in range(3)] for j in range(len(atoms))]
        atoms.translate(tr2d)

    elif dim == 1:
        pass
    elif dim == 2:
        if structure == None:
            print(f"input structure name using -n ")
            sys.exit(1)
        elif sturucture == 'graphene':
            atoms = bd.graphene(size=(4,4,1), vacuum=vac)       # does not make z-axis
            view(image)
        elif structure == 'slab':
            pass
    elif dim == 'bulk':
        if name in metals.keys():
            if metals[name] == 'fcc':
                from ase.lattice.cubic import FaceCenteredCubic
                atoms = FaceCenteredCubic(symbol=name, size=size, pbc=True)
        elif status == 'gas':
            #coord = make_random_coordinates(natom, latt)
            coord = make_ordered_coordinates(natom, latt)
            #print(coord)
            atoms = Atoms(symbols='Ar256', positions=coord, cell=(latt,latt,latt), pbc=True)
        
    return atoms

def main():
    parser = argparse.ArgumentParser('ASE builder for structure')
    #parser.add_argument('module', default='build', choices=['build'], help='ASE jobs')
    parser.add_argument('-suf', '--out_suffix', help='POSCAR.suff for output filename')
    parser.add_argument('-d', '--dim', default=2, type=int, choices=[0,1,2,3], help='basic structure of system')
    parser.add_argument('-s', '--structure', help='structure name following dimension')
    parser.add_argument('-n', '--name', help='structure name')
    parser.add_argument('-si', '--size', nargs='*', default=[1], type=int, help='size of supercell')
    parser.add_argument('-v', '--vacuum', nargs='*', default=[10.0], type=float, help='size of supercell')
    parser.add_argument('-ak', '--atomkind', default='Pt', help='atom species')
    parser.add_argument('-u', '--usage', action='store_true', help='shows usage and exit')
    args = parser.parse_args()

    Usage.cluster="ase_build.py -d 0 -s icosa -ak Pd -v 20 -si 2 -suf test : Pd icosahedron cluster of size 2\
                    \nase_build.py -d 0 -s deca  -ak Pd -v 20 -si 2 2\t\t: Pd decahedron with length 2 & height 2\
                    \nase_build.py -d 0 -s octa  -ak Pd -v 20 -si 3  \t\t: Pd octahedron with size 3 -> truncated octa"

    if args.usage:
        print(Usage.cluster)
    else:
        image = build_structure(args.dim, args.structure, args.name, args.size, args.vacuum, args.atomkind)
        ### make output filename with suffix after POSCAR
        suff=''
        if args.atomkind:
            suff += args.atomkind
        if args.structure:
            suff += args.structure
        if args.size:
            st = list2str(args.size)
            suff += st
        ### if suffix exists, use it
        if args.out_suffix:
            outf = 'POSCAR.'+args.out_suffix
        elif suff:
            outf = 'POSCAR.'+suff
        else:
            outf = 'POSCAR'
        print(f"suffix {suff}: write the image to {outf}")
        write(outf, image, format='vasp')

    return 0

if __name__ == '__main__':
    main()
