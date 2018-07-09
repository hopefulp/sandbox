#!/usr/bin/python

import argparse
import os

def main():

    parser = argparse.ArgumentParser(description='prepare vasp input files')
    parser.add_argument('-o', '--odir', help='copy from old directory')
    parser.add_argument('new_dir', help='make new directory')
    parser.add_argument('-f', '--files', nargs='*', help='add files to be copied')

    args = parser.parse_args()

    files = ['INCAR', 'POSCAR', 'POTCAR', 'KPOINTS']

    if args.files:
        files.extend(args.files)

    cmd = 'mkdir ' + args.new_dir
    os.system(cmd)

    for file in files:
        if args.odir:
            ofile = args.odir + '/' + file
        else:
            ofile = file
        cmd =   'cp ' + ofile + ' ' + args.new_dir
        os.system(cmd)
    return 0

if __name__ == '__main__':
    main()
