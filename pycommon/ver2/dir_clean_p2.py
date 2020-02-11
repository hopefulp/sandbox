#!/usr/bin/python

import argparse
import os
import re
import sys
from common_p2 import *

q_list=[]

def d_clean(dirs,work,prefix, suffix, matches, exclude,excl_fnames, linux_job,new_dir,Lshowmatch):

    pwd = os.getcwd()
    if len(dirs) > 1:
        print "use only one directory"
        sys.exit(1)
    else:
        d = pwd + '/' + dirs[0]

    matches=[]
    if work == None:
        if prefix:
            f_list=get_files_prefix(prefix, d)
        if suffix:
            f_list=get_files_suffix(suffix, d)
        if matches:
            f_list=get_files_match(matches, d,Lshowmatch)
        if exclude:
            f_list=get_files_exclude(exclude, d)
        if f_list and excl_fnames:
            for f in excl_fnames:
                if f in f_list:
                    f_list.remove(f)
    elif work == 'qchem':
        pass
    elif work == 'vasp':
        f_list=os.listdir(d)
        if excl_fnames:
            for f in f_list:
                if os.access(d+'/'+f, os.X_OK):
                    excl_fnames.append(f)
            for efile in excl_fnames:
                if efile in f_list:
                    f_list.remove(efile)

    elif work == 'pbs':
        matches=['\.e\d', '\.o\d', '\.pe\d', '\.po\d']
        f_list = get_files_match(matches, d, Lshowmatch)
            

    f_list.sort()
    ### Make command list
    for f in f_list:
        fname = d+'/'+f
        comm = "%s %s" % (linux_job, fname)
        if linux_job == 'mv':
            comm += " %s" % new_dir
        print comm
        q_list.append(comm)
        
    #print "all %s files" % len(f_list)
    ### Show command list and run
    if f_list:
        q = "will you %s %s files? " % (linux_job, len(f_list))
        if yes_or_no(q):
            i = 0
            for comm in q_list:
                os.system(comm)
                i += 1
            print "%s files are removed" % i         

    return 0


def main():
    parser = argparse.ArgumentParser(description='to clean directory in qchem')
    parser.add_argument('dirs', default='.', nargs='+', help='input work directories')
    parser.add_argument('-w', '--work', choices=['qchem','ai','vasp','pbs'],help='remove depending on job')
    parser.add_argument('-p', '--prefix', nargs='*', help='remove with prefix')
    parser.add_argument('-s', '--suffix', nargs='*', help='remove with suffix')
    parser.add_argument('-m', '--match', nargs='*', help='remove matching file')
    parser.add_argument('-e', '--exclude', nargs='*', help='remove all files except list') 
    parser.add_argument('-ef', '--excluded_files', nargs='*', help='save this file') 
    parser.add_argument('-j', '--job', default='rm', choices=['rm','mv'], help='how to treat files')
    parser.add_argument('-jd', '--mv_dir', default='tmp', help='directory where files to move')
    parser.add_argument('-ms', '--match_show', action='store_true')
    args = parser.parse_args()

    if args.work==None and args.prefix==None and args.suffix==None and args.match==None:
        print "input -w|-p|-s|-m"
        print "use %s -h for help" % os.path.basename(__file__)
        sys.exit(0)
    if args.work == 'vasp' and not args.excluded_files:
        args.excluded_files=['POSCAR','POTCAR','KPOINTS','INCAR']
            
    d_clean(args.dirs,args.work,args.prefix,args.suffix,args.match,args.exclude,args.excluded_files,args.job,args.mv_dir,args.match_show)
    return 0

if __name__ == '__main__':
    main()
