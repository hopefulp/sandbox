#!/usr/bin/python

import argparse
import os
import re
from common import *

q_list=[]

def clean(prefix, suffix, matches, exclude):

    pwd = os.getcwd()
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
    parser.add_argument('-p', '--prefix', nargs='*', help='remove with prefix')
    parser.add_argument('-s', '--suffix', nargs='*', help='remove with suffix')
    parser.add_argument('-m', '--match', nargs='*', help='remove matching file')
    #parser.add_argument('-r', '--reverse', nargs='*', help='remove matching file')
    #parser.add_argument('-e', '--exclude', default=["out","rem"], nargs='+', help='remove all files except list') 
    parser.add_argument('-e', '--exclude', nargs='*', help='remove all files except list') 
    args = parser.parse_args()

    clean(args.prefix, args.suffix, args.match, args.exclude)
    return 0

if __name__ == '__main__':
    main()
