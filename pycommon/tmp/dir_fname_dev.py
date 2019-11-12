#!/usr/bin/python

import argparse
import re
import os
import sys
from common_p2 import *

def cmd(job, m_tag, pattern,L_dir, exceptions, dir_destination,rn_type,new_name,run):
    pwd = os.getcwd()
    # pattern should have 1 element
    print pattern
    l_file = get_files_pattern(m_tag, pattern, pwd,L_dir)
    n=0
    #print("\n".join(l_file))                       
    if exceptions:
        for fex in exceptions:
            l_file.remove(fex)
    # run job for the file list
    #print("\n".join(l_file))
    print len(l_file), " files are selected"
    if job == "ls":
        if l_file:
            print("\n".join(l_file))
    elif job == "mvdir" :
        if not os.path.isdir(dir_destination):
            comm = "mkdir " + dir_destination
            os.system(comm)
        for file in l_file:
            comm = job + " " + file + " " + dir_destination
            if run:
                os.system(comm)
            else:
                print comm
    elif job == 'rename':
        if not new_name: 
            print "add new name with -a or "
            exit(10)
        else:
            for fname in l_file:
                #new_name = fname.replace(pattern[0], new_pattern)
                if re.search('\.', fname):
                    f_name = re.split('\.', fname)
                #print f_name
                    new_file=f_name[0]+new_name+'.'+f_name[1]
                else:
                    new_file=fname+new_name
                comm  = "mv  " + fname + " " + new_file
                if run:
                    os.system(comm)
                else:
                    print comm
    elif job == 'rm':
        for file in l_file:
            comm = job + " " + file
            if run :
                os.system(comm)
            else:
                print comm
    
    if(job == 'rm' or job == 'mvdir' or job == 'rename') and run == False:
        print "use '-r' for run"
    return                    


def main():
    parser = argparse.ArgumentParser(description='directory management')
    parser.add_argument( 'job', choices=['ls', 'mvdir', 'rm', 'rename'],  help='shell command')
    group = parser.add_mutually_exclusive_group()
    group.add_argument( '-p', '--prefix', nargs='*', help='prefix of filename')
    group.add_argument( '-s', '--suffix', nargs='*', help='list several suffixes')
    group.add_argument( '-m', '--match', nargs='*', help='find matching string')
    parser.add_argument( '-i', '--include_dir', action='store_true', help='include dirname to filename')
    parser.add_argument( '-xcl', '--excluded', nargs='*', help='filename to be excluded')
    parser.add_argument( '-d', '--directory', type=str, default='tmppy', help='target directory to move files')
    #parser.add_argument( '-rn', '--rename', help='rename filename')
    group_rn=parser.add_mutually_exclusive_group()
    group_rn.add_argument('-a', '--append', help="add suffix by -a to the original filename without extension")
    parser.add_argument( '-r', '--run', action='store_true', help='run or not-False')
    args = parser.parse_args()

    if args.prefix:
        matching=args.prefix
        m_tag = 'p'
    elif args.suffix:
        matching=args.suffix
        m_tag = 's'
    elif args.match:
        matching=args.match
        m_tag = 'm'
    else:
        print "matching should be given"
        return 1

    if args.job == "rename":
        if args.append:
            new_name = args.append
            rn_type = 'a'   # for append except extension
    #cmd(args.job, m_tag, matching, args.excluded, args.directory,args.rename,rn_type,new_name,args.run)
    cmd(args.job, m_tag, matching,args.include_dir, args.excluded, args.directory,rn_type,new_name,args.run)

if __name__ == "__main__":
    main()
