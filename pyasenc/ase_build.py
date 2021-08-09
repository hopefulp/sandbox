#!/home/joonho/anaconda3/bin/python

import argparse
from ase import Atoms, Atom, units
import ase.build as bd
from ase.visualize import view
from ase.io.formats import read, iread, write, string2index
from random import uniform
import numpy as np

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

    for i in range(nx):
            



def build_structure(dim, status=None, name=None, natom=10, latt=50):
    if dim == 'mol':
        pass
    elif dim == '1d':
        pass
    elif dim == '2d':
        if name == 'graphene':
            image = bd.graphene(size=(4,4,1), vacuum=10.0)       # does not make z-axis
            view(image)
    elif dim == 'slab':
        pass
    elif dim == 'bulk':
        if status == 'gas':
            #coord = make_random_coordinates(natom, latt)
            coord = make_ordered_coordinates(natom, latt)
            #print(coord)
            image = Atoms(symbols='Ar256', positions=coord, cell=(latt,latt,latt), pbc=True)
        
    return image 

def main():
    parser = argparse.ArgumentParser('ASE builder for structure')
    parser.add_argument('module', default='build', choices=['build'], help='ASE jobs')
    parser.add_argument('-d', '--dim', default='2d', choices=['mol', '1d', '2d', 'slab', 'bulk'], help='basic structure of system')
    parser.add_argument('-n', '--name', help='structure name')
    args = parser.parse_args()

    if args.module == 'build':
        image = build_structure(args.dim, args.name)
        write('POSCAR', image, format='vasp')

    return 0

if __name__ == '__main__':
    main()
