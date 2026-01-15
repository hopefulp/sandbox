#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
import sys
from common import *
from libdir import walk_dirs,  dir_clean, dir_touch, CleanConfig
from functools import partial

def main():
    parser = argparse.ArgumentParser( description='directory operations (clean, touch, ...)' )
    # NEW: top-level job
    parser.add_argument('-j', '--job', choices=['clean', 'touch'], default='clean', help='directory job type' )
    # NEW: subjob (old -j)
    parser.add_argument( '-sj', '--shelljob', choices=['rm', 'mv', 'cp', 'ln'], default='rm', help='file operation for clean' )
    parser.add_argument('-d', '--dir1', help='input work directory')
    ### choose: p, s, m, w
    gselection = parser.add_mutually_exclusive_group()
    gselection.add_argument( '-w', '--works', nargs='*', choices=['qchem','amp','vasp','pbs','slurm','lmp','nc'], help='remove depending on job')
    gselection.add_argument('-p', '--prefix', nargs='*')
    gselection.add_argument('-s', '--suffix', nargs='*')
    gselection.add_argument('-m', '--match', nargs='*')
    parser.add_argument('-e', '--exclude', nargs='*')
    parser.add_argument('-ef', '--excluded_files', nargs='*')       # for filename to be excluded
    parser.add_argument('-jd', '--new_dir', default='tmp')
    parser.add_argument('-ms', '--match_show', action='store_true')
    parser.add_argument('-a', '--all_remove', action='store_true')
    parser.add_argument('-y', '--yes', action='store_true')

    args = parser.parse_args()

    # ---- validation ----
    if args.job == 'clean':
        if args.works:
            selection = 'w'
            slist = args.works
        elif args.prefix:
            selection = 'p'
            slist = args.prefix
        elif args.suffix:
            selection = 's'
            slist = args.suffix
        elif args.match:
            selection = 'm'
            slist = args.match
        else:
            parser.error("one of -w/-p/-m/-s is required for clean")

        if args.works and 'vasp' in args.works and not args.excluded_files:
            args.excluded_files = ['POSCAR', 'POTCAR', 'KPOINTS', 'INCAR']

        config = CleanConfig(
            selection=selection,
            slist=slist,
            works=args.works,
            shelljob=args.shelljob,
            exclude=args.exclude,
            excl_fnames=args.excluded_files,
            new_dir=args.new_dir,
            show_match=args.match_show,
            all_rm=args.all_remove,
            yes=args.yes,
            )

        clean_work = partial(dir_clean, config=config)
        walk_dirs(args.dir1, clean_work)
    elif args.job == 'touch':
        walk_dirs(args.dir1, dir_touch)
    else:
        print("Present dir_works: clean, touch inserted by -j")
        sys.exit(0)
    return 0

if __name__ == '__main__':
    main()
