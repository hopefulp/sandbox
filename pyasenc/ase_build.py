#!/home/joonho/anaconda3/bin/python

import argparse
import ase.build as bd
from ase.visualize import view
from ase.io.formats import read, iread, write, string2index

def build_structure(dim, name):
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
        pass
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
