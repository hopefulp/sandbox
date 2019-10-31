#!/usr/bin/python

import argparse
import os
import re
import sys
from common_p2 import *

q_list=[]

def d_clean(work,w_option,prefix, suffix, matches, exclude,excl_fnames, linux_job,new_dir):

    pwd = os.getcwd()
    if work == None:
        if prefix:
            f_list=get_files_prefix(prefix, pwd)
        if suffix:
            f_list=get_files_suffix(suffix, pwd)
        if matches:
            f_list=get_files_match(matches, pwd)
        if exclude:
            f_list=get_files_exclude(exclude, pwd)
        if f_list and excl_fnames:
            for f in excl_fnames:
                if f in f_list:
                    f_list.remove(f)
    elif work == 'qchem':
        if w_option == None:
            print "work==qchem requires pid number with -wo n"
            sys.exit(2)
        matches=[]
        matches.append(".e"+str(w_option))
        matches.append(".pe"+str(w_option))
        matches.append(".o"+str(w_option))
        matches.append(".pp"+str(w_option))
        f_list=get_files_match(matches, pwd)
    
    elif work == 'vasp':
        f_list=os.listdir(pwd)
        if excl_fnames:
            for efile in excl_fnames:
                if efile in f_list:
                    f_list.remove(efile)

    f_list.sort()
    for f in f_list:
        comm = "%s %s" % (linux_job, f)
        if linux_job == 'mv':
            comm += " %s" % new_dir
        print comm
        q_list.append(comm)
        
    #print "all %s files" % len(f_list)
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
    parser.add_argument('-w', '--work', choices=['qchem','ai','vasp'],help='remove depending on job')
    parser.add_argument('-wo', '--work_option', type=int, help='Q-Chem: parallel output number')
    parser.add_argument('-p', '--prefix', nargs='*', help='remove with prefix')
    parser.add_argument('-s', '--suffix', nargs='*', help='remove with suffix')
    parser.add_argument('-m', '--match', nargs='*', help='remove matching file')
    parser.add_argument('-e', '--exclude', nargs='*', help='remove all files except list') 
    parser.add_argument('-ef', '--excluded_files', nargs='*', help='save this file') 
    parser.add_argument('-j', '--job', default='rm', choices=['rm','mv'], help='how to treat files')
    parser.add_argument('-jd', '--mv_dir', default='tmp', help='directory where files to move')
    args = parser.parse_args()

    if args.work==None and args.prefix==None and args.suffix==None and args.match==None:
        print "input -w|-p|-s|-m"
        print "use %s -h for help" % os.path.basename(__file__)
        sys.exit(0)
    if args.work == 'vasp' and not args.excluded_files:
        args.excluded_files=['POSCAR','POTCAR','KPOINTS','INCAR']
            
    d_clean(args.work,args.work_option,args.prefix,args.suffix,args.match,args.exclude,args.excluded_files,args.job,args.mv_dir)
    return 0

if __name__ == '__main__':
    main()
