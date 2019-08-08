#!/usr/bin/python

import argparse
import os
import re
import sys
from common_p2 import *

def search_dir(job, dir, foption, lfile_exc, exe_tag):
    os.chdir(dir)
    print("####.... enter %s directory" % dir)
    p_files = os.listdir('.')
    #print("%s" % type(foption))
    for file in p_files:
        if not os.path.isdir(file):
            if file in lfile_exc:
                pass
            else:
                if isinstance(foption, list):
                    if file in foption:
                        cmd = job + ' ' + file 
                        print cmd
                        if exe_tag:
                            os.system(cmd)
                            print("%s was removed" % file)
                    else:
                        pass
                elif isinstance(foption, float):
                    fsize = os.path.getsize(file)
                    fsize /= 1024.0 * 1024.0
                    #print("%f %f" % (foption, fsize))
                    if fsize > foption:
                        cmd = job + ' ' + file
                        print cmd
                        if exe_tag:
                            os.system(cmd)
                            print("%s was removed" % file)
                    else:
                        pass
        # if directory
        else:
            search_dir(job, file, foption, lfile_exc, exe_tag)
    os.chdir('..')                
    print "####....... exit %s directory" % dir

def f_clean(dirs,work,work_suffix,linux_job,mv_dir,files,option,opt_str,exclude,excluded_files):
    if not dirs:
        pwd = os.getcwd()
    if work == None:
        if option == 'p':
            get_files = get_files_prefix
        elif option == 's':
            get_files = get_files_suffix
        elif option == 'm':
            get_files = get_files_match
        else:
            print "error in option and work"
        

            #f_list=get_files_exclude(exclude, pwd)
        #if f_list and excl_fnames:
        #    for f in excl_fnames:
        #        if f in f_list:
        #            f_list.remove(f)
    elif work == 'qchem':
        default_str=['\.e', '\.pe', '\.o', '\.po']
        opt_str=[]
        if not work_suffix:
            print "input work suffix number for %s with -ws" % work
            sys.exit(10)
        else:
            for st in default_str:
                st+=work_suffix
                opt_str.append(st)

        get_files = get_files_match
    else:
        print "No %s" % work
        sys.exit(1)

    q_list=[]
    f_list=get_files(opt_str, pwd)
    #print f_list
    #sys.exit(5)

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

    parser = argparse.ArgumentParser(description='remove files recursively u. size or filename')
    parser.add_argument('-w', '--work', choices=['qchem','ai'],help='remove depending on job')
    parser.add_argument('-ws', '--work_suf', help='extra keywords for work')
    parser.add_argument('-j', '--job', default='rm', choices=['rm','mv'], help='how to treat files')
    parser.add_argument('-jd', '--mv_dir', default='tmp', help='directory where files to move')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-f', '--files', nargs='*', help='files to be remode')
    group.add_argument('-o', '--option', default=['p','s','m','z'], help='p for prefix, s for suffix, m for middle, z for size')

    parser.add_argument('-of', '--opt_str', help='common key word for files/minimum size of files')
    parser.add_argument('-e', '--exclude', nargs='*', help='remove all files except list')
    parser.add_argument('-ef', '--excluded_files', nargs='*', help='save this file')

    parser.add_argument('-d', '--dirs', nargs='*', help='select directorys')
    parser.add_argument('-r', '--recursive', action='store_true', help='search the sub-directories')
    args = parser.parse_args()

    job=args.job

    #parser.add_argument('-x', '--xfiles', nargs='*', default=['INCAR', 'POSCAR', 'POTCAR', 'KPOINTS', 'IBZKPT', 'CONTCAR'], help='excluded files')
    del_files = ['WAVECAR', 'CHGCAR', 'CHG']
    #print args

    if args.recursive: 
        if args.files:
            foption = args.files
        else:
            foption = args.size

        for dir in args.dir:
            search_dir(job, dirs, foption, args.xfiles, args.exe)
    else:
        if args.work==None and args.option==None:
            print "input -w|-p|-s|-m"
            print "use %s -h for help" % os.path.basename(__file__)
            sys.exit(0)
        else:
            f_clean(args.dirs,args.work,args.work_suf,args.job,args.mv_dir,args.files,args.option,args.opt_str,args.exclude,args.excluded_files)

    return 0

if __name__ == '__main__':
    main()
