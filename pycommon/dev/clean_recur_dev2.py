#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
import sys
from common import *
from clean1d import dir_clean

def d_clean_recursive(dir1,works,linux_job,prefix, suffix, matches, exclude,excl_fnames,new_dir,Lshowmatch,Lall_rm, Lyes):
    print(f"{dir1}")
    if not dir1:
        pwd = os.getcwd()
        dir1 = pwd
    print(f"####.... enter {dir1} directory")
    os.chdir(dir1)
    ### call d_clean in clean1d.py
    dir_clean(None,works,linux_job,prefix, suffix, matches, exclude,excl_fnames,new_dir,Lshowmatch,Lall_rm, Lyes)
    p_files = os.listdir('.')
    for f in p_files:
        ### without check not link-directory, this makes error
        if os.path.isdir(f) and not os.path.islink(f):
            d_clean_recursive(f,works,linux_job,prefix, suffix, matches, exclude,excl_fnames,new_dir,Lshowmatch,Lall_rm, Lyes)

    os.chdir('..')
    print(f"####....... exit {dir1} directory")
    return 0


def main():
    parser = argparse.ArgumentParser(description='to clean directory in qchem')
    parser.add_argument('-d', '--dir1',  help='input work directories')
    parser.add_argument('-w', '--works', nargs='+', choices=['qchem','amp','vasp','pbs','lmp'],help='remove depending on job')
    parser.add_argument('-j', '--job', default='rm', choices=['rm','mv','cp','ln'], help='how to treat files [rm|cp|mv]')
    parser.add_argument('-p', '--prefix', nargs='*', help='remove with prefix')
    parser.add_argument('-s', '--suffix', nargs='*', help='remove with suffix')
    parser.add_argument('-m', '--match', nargs='*', help='remove matching file')
    parser.add_argument('-e', '--exclude', nargs='*', help='remove all files except list') 
    parser.add_argument('-ef', '--excluded_files', nargs='*', help='save this file') 
    parser.add_argument('-jd', '--new_dir', default='tmp', help='directory where files to move')
    parser.add_argument('-ms', '--match_show', action='store_true')
    parser.add_argument('-a', '--all_remove', action='store_true', help='remove all the files')
    parser.add_argument('-y', '--yes', action='store_true', help='execute command')
    args = parser.parse_args()

    if args.works==None and args.prefix==None and args.suffix==None and args.match==None:
        print("input -w|-p|-s|-m")
        print("use %s -h for help" % os.path.basename(__file__))
        sys.exit(0)
    if 'vasp' in args.works and not args.excluded_files:
        args.excluded_files=['POSCAR','POTCAR','KPOINTS','INCAR']
    #if args.work == 'amp' and not args.excluded_files:
    #    args.excluded
    d_clean_recursive(args.dir1,args.works,args.job,args.prefix,args.suffix,args.match,args.exclude,args.excluded_files,args.new_dir,args.match_show,args.all_remove, args.yes)
    return 0

if __name__ == '__main__':
    main()
