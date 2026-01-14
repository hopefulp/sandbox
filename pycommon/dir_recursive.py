#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
import sys
from common import *
from libdir import walk_dirs,  dir_clean, dir_touch
from functools import partial

def main():
    parser = argparse.ArgumentParser( description='directory operations (clean, touch, ...)' )
    # NEW: top-level job
    parser.add_argument('-j', '--job', choices=['clean', 'touch'], default='clean', help='directory job type' )
    # NEW: subjob (old -j)
    parser.add_argument( '-sj', '--subjob', choices=['rm', 'mv', 'cp', 'ln'], default='rm', help='file operation for clean' )
    parser.add_argument('-d', '--dir1', help='input work directory')
    parser.add_argument( '-w', '--works', nargs='+', choices=['qchem', 'amp', 'vasp', 'pbs', 'lmp'], help='remove depending on job')
    parser.add_argument('-p', '--prefix', nargs='*')
    parser.add_argument('-s', '--suffix', nargs='*')
    parser.add_argument('-m', '--match', nargs='*')
    parser.add_argument('-e', '--exclude', nargs='*')
    parser.add_argument('-ef', '--excluded_files', nargs='*')
    parser.add_argument('-jd', '--new_dir', default='tmp')
    parser.add_argument('-ms', '--match_show', action='store_true')
    parser.add_argument('-a', '--all_remove', action='store_true')
    parser.add_argument('-y', '--yes', action='store_true')

    args = parser.parse_args()


    # ---- validation ----
    if args.job == 'clean':
        if (args.works is None and args.prefix is None and args.suffix is None and args.match is None):
            print("input -w|-p|-s|-m")
            print(f"use {os.path.basename(__file__)} -h for help")
            sys.exit(10)
        if args.works and 'vasp' in args.works and not args.excluded_files:
            args.excluded_files = ['POSCAR', 'POTCAR', 'KPOINTS', 'INCAR']
        
        clean_work = partial( dir_clean, args.works, args.subjob,      # <-- old job goes here
            args.prefix, args.suffix, args.match, args.exclude, args.excluded_files, args.new_dir, args.match_show,
            args.all_remove, args.yes
        )
    #d_clean_recursive(args.dir1,args.works,args.job,args.prefix,args.suffix,args.match,args.exclude,args.excluded_files,args.new_dir,args.match_show,args.all_remove, args.yes)
        walk_dirs(args.dir1, clean_work)
    elif args.job == 'touch':
        walk_dirs(args.dir1, dir_touch)
        # later:
        # d_touch_recursive(args.dir1, ...)
    else:
        print("Present dir_works: clean, touch inserted by -j")
        sys.exit(0)
    return 0

if __name__ == '__main__':
    main()
