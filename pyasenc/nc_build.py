#!/home/joonho/anaconda3/bin/python

import argparse
import nanocore as nc
from nanocore import io, carbonlab
from nanocore.vasp import Vasp

def print_image(image):
    pos = image.get_positions()
    cell = image.get_cell()
    natom = len(pos)
    print(f"cell {cell}")
    print(f"natoms {natom}")

    print(*pos, sep='\n')
    return 0

def build_structure(dim, name, subname, sc, dz, vac, Lvis, Lsave, ofname, fformat):
    if dim == 'mol':
        pass
    elif dim == '1d':
        pass
    elif dim == '2d':
        if name == 'graphene':
            if subname == 'hexa':
                image = carbonlab.grp(sc[0], sc[1], vacuum=vac)       # does not make z-axis
            elif subname == 'rect':
                image = carbonlab.grp_rect(sc[0], sc[1], vacuum=vac)       # does not make z-axis
            elif subname == 'zz':
                image = carbonlab.gnr_zz(sc[0], sc[1], vacuum=vac)       # does not make z-axis
            elif subname == 'ac':
                image = carbonlab.gnr_ac(sc[0], sc[1], vacuum=vac)       # does not make z-axis
        image.select_all()
        trans = [0, 0, dz*vac]
        image.translate(*trans)
    elif dim == 'slab':
        pass
    elif dim == 'bulk':
        pass
    if Lvis:
        from nanocore import vis
        vis.show_xcrysden(image)
    if Lsave:
        if fformat == 'vasp':
            print("make Vasp object")
            vas = Vasp(image)
            fname = 'POSCAR'
            if ofname:
                fname += '.'+ofname
            vas.write_POSCAR(file_name=fname)
            print(f"wrote poscar to {fname}")
            
    return 0 

def main():
    parser = argparse.ArgumentParser('NanoCore builder for structure')
    parser.add_argument('module', default='build', choices=['build'], help='NonoCore jobs')
    parser.add_argument('-d', '--dim', default='2d', choices=['mol', '1d', '2d', 'slab', 'bulk'], help='basic structure of system')
    parser.add_argument('-n', '--name', help='structure name')
    parser.add_argument('-sn', '--subname', default='hexa', help='structure subname')
    parser.add_argument('-sc', '--supercell', nargs='+', default=[3,4,1], type=int, help='size of supercell')
    parser.add_argument('-t', '--translate', default=0.4, help='shift of plane in z-axis')
    parser.add_argument('-v', '--vacuum', default=10.0, help='shift of plane in z-axis')
    parser.add_argument('-g', '--visual', action='store_true', help='view structure via xcrysden')
    parser.add_argument('-s', '--save', action='store_true', help='write built structure to a file format')
    parser.add_argument('-o', '--outfile', default='tmp', help='output filename')
    parser.add_argument('-fm', '--fformat', default='vasp', help='output file format')

    args = parser.parse_args()

    if args.module == 'build':
        build_structure(args.dim, args.name, args.subname, args.supercell, args.translate, args.vacuum, args.visual, args.save, args.outfile, args.fformat)

    return 0

if __name__ == '__main__':
    main()
