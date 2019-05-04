#!/usr/bin/python

import argparse
import os
import re
import sys
from common_p2 import *

q_list=[]

def f_clean(job,subjob,prefix, suffix, matches, exclude):

    pwd = os.getcwd()
    if job == None:
        if prefix:
            f_list=get_files_prefix(prefix, pwd)
            for file in f_list:
                comm = "rm %s" % file
                print comm
                q_list.append(comm)
        if suffix:
            f_list=get_files_suffix(suffix, pwd)
            for file in f_list:
                comm = "rm %s" % file
                print comm
                q_list.append(comm)
        if matches:
            f_list=get_files_match(matches, pwd)
            for file in f_list:
                comm = "rm %s" % file
                print comm
                q_list.append(comm)
        if exclude:
            f_list=get_files_exclude(exclude, pwd)
            for file in f_list:
                comm = "rm %s" % file
                print comm
                q_list.append(comm)
    elif job == 'qchem':
        if subjob == None:
            print "job==qchem requires pid number with -js n"
            sys.exit(2)
        matches=[]
        matches.append(".e"+str(subjob))
        matches.append(".pe"+str(subjob))
        matches.append(".o"+str(subjob))
        matches.append(".pp"+str(subjob))

        f_list=get_files_match(matches, pwd)
        for file in f_list:
            comm = "rm %s" % file
            print comm
            q_list.append(comm)
        

    #print "all %s files" % len(f_list)            
    q = "will you remove %s files? " % len(f_list)
    if yes_or_no(q):
        i = 0
        for comm in q_list:
            os.system(comm)
            i += 1
        print "%s files are removed" % i         

    return 0


def main():
    parser = argparse.ArgumentParser(description='to clean directory in qchem')
    parser.add_argument('-j', '--job', choices=['qchem','ai'],help='remove depending on job')
    parser.add_argument('-js', '--sub_job', type=int, help='Q-Chem: parallel output number')
    parser.add_argument('-p', '--prefix', nargs='*', help='remove with prefix')
    parser.add_argument('-s', '--suffix', nargs='*', help='remove with suffix')
    parser.add_argument('-m', '--match', nargs='*', help='remove matching file')
    parser.add_argument('-e', '--exclude', nargs='*', help='remove all files except list') 
    args = parser.parse_args()

    if args.job==None and args.prefix==None and args.suffix==None and args.match==None:
        print "input -j|-p|-s|-m"
        print "use %s -h for help" % os.path.basename(__file__)
        sys.exit(0)

    f_clean(args.job,args.sub_job,args.prefix,args.suffix,args.match,args.exclude)
    return 0

if __name__ == '__main__':
    main()
