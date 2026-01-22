#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
import sys
from common import *
from libdir.walker import walk_dirs
from libdir.clean import dir_clean, CleanConfig
from libdir.jobs import dir_touch
from functools import partial

def main():
    parser = argparse.ArgumentParser( description='directory operations (clean, touch, ...)' )
    # NEW: top-level dirjob
    parser.add_argument('-j', '--dirjob', choices=['clean', 'touch'], default='clean', help='directory dirjob type' )
    # NEW: subjob (old -j)
    parser.add_argument( '-sj', '--shelljob', choices=['rm', 'mv', 'cp', 'ln'], default='rm', help='file operation for clean' )
    parser.add_argument('-d', '--dir1', help='input work directory')
    ### choose: p, s, m, w
    gselection = parser.add_mutually_exclusive_group()
    gselection.add_argument( '-w', '--works', nargs='*', choices=['qchem','amp','vasp','pbs','slurm','lmp','nc'], help='remove depending on dirjob')
    gselection.add_argument('-p', '--prefix', nargs='*')
    gselection.add_argument('-s', '--suffix', nargs='*')
    gselection.add_argument('-m', '--match', nargs='*')
    parser.add_argument('-e', '--exclude', nargs='*')
    parser.add_argument('-ef', '--excluded_files', nargs='*')       # for filename to be excluded
    parser.add_argument('-jd', '--new_dir', default='tmp')
    parser.add_argument('-ms', '--match_show', action='store_true')
    parser.add_argument('-a', '--all_remove', action='store_true')
    parser.add_argument('-y', '--yes', action='store_true')
    parser.add_argument('-R', '--recursive', action='store_true')
    parser.add_argument('-u', '--usage', action='store_true')

    args = parser.parse_args()

    if args.usage:
        print(f"1. select one of jobs [clean|touch")
        print(f"2. select one of selection rule [-p|-m|-s|-w]")
        print(f"3. -R for recursive")
        print(f"Usage::")
        print(f"\tdir_jobs.py -j clean -w pbs")
        print(f"\tdir_jobs.py -j touch -R")
        sys.exit(0)
    path = os.path.abspath(args.dir1 or os.getcwd())

    if args.dirjob == 'clean':

        if args.prefix:
            selection = 'p'
            slist = args.prefix
        elif args.suffix:
            selection = 's'
            slist = args.suffix
        elif args.match:
            selection = 'm'
            slist = args.match
        elif args.works:
            selection = 'w'
            slist = args.works        
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

        if args.recursive:
            walk_dirs(args.dir1, clean_work)
        else:
            clean_work(path)
    elif args.dirjob == 'touch':
        if args.recursive:
            walk_dirs(args.dir1, dir_touch)
        else:
            dir_touch(path)
    else:
        print("Present dir_works: clean, touch inserted by -j")
        sys.exit(0)
    return 0

if __name__ == '__main__':
    main()
