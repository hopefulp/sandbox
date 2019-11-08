#!/usr/bin/python

import argparse
import re
import os
import sys
from common import *

def cmd(job, m_tag, pattern, exceptions, dir_out,rn_type,new_name,run):
    pwd = os.getcwd()
    # pattern should have 1 element
    print pattern
    if m_tag != 'f':
        l_file = get_files_pattern(m_tag, pattern, pwd)
    else:
        l_file = pattern
    n=0

    #print("\n".join(l_file))                       
    if exceptions:
        for fex in exceptions:
            l_file.remove(fex)
    # run job for the file list
    #print("\n".join(l_file))
    print len(l_file), " files are selected"
    command=[]
    if job == "ls":
        if l_file:
            print("\n".join(l_file))
    elif job == "mvdir" :
        if not os.path.isdir(dir_out):
            comm = "mkdir " + dir_out
            os.system(comm)
        for f in l_file:
            comm = job + " " + f + " " + dir_out
            print comm
    elif job == 'rename':
        if not new_name: 
            print "add new name with -a or "
            exit(10)
        else:
            for fname in l_file:
                #new_name = fname.replace(pattern[0], new_pattern)
                f_name = re.split('\.', fname)
                #print f_name
                new_file=f_name[0]+new_name+'.'+f_name[1]
                comm  = "mv  " + fname + " " + new_file
                if run:
                    os.system(comm)
                else:
                    print comm
    elif job == 'rm':
        for f in l_file:
            comm = job + " " + f
            if run :
                os.system(comm)
            else:
                print comm
    elif job == 'cp':
        
        for f in l_file:
            for d in dir_out:
                com = 'cp %s %d' % (f, d)
                commands.append(com)

                
    
    if(job == 'rm' or job == 'mvdir' or job == 'rename') and run == False:
        print "use '-r' for run"
    for com in commands:
        print com
    q = "will you run? "
    if yes_or_no(q):
        for com in commands:
            os.system(comxm)
    return                    


def main():
    parser = argparse.ArgumentParser(description='Command Line Interface to deal with directory')
    parser.add_argument( 'job', choices=['ls', 'mvdir', 'rm', 'rename', 'cp'],  help='shell command')
    group = parser.add_mutually_exclusive_group()
    group.add_argument( '-p', '--prefix', nargs='*', help='prefix of filename')
    group.add_argument( '-s', '--suffix', nargs='*', help='list several suffixes')
    group.add_argument( '-m', '--match', nargs='*', help='find matching string')
    group.add_argument( '-f', '--files', nargs='*', help='input file list')
    parser.add_argument( '-xcl', '--excluded', nargs='*', help='filename to be excluded')
    parser.add_argument( '-d', '--directory', args='+', type=str, default='tmppy', help='target directory to move files')
    #parser.add_argument( '-rn', '--rename', help='rename filename')
    group_rn=parser.add_mutually_exclusive_group()
    group_rn.add_argument('-a', '--append', default='_new', help="add '_new' to the original filename withoug extension")
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
    elif args.ifile:
        matching=args.files
        m_tag = 'f'
    else:
        print "matching should be given"
        return 1

    if args.job == "rename":
        if args.append:
            new_name = args.append
            rn_type = 'a'   # for append except extension
    #cmd(args.job, m_tag, matching, args.excluded, args.directory,args.rename,rn_type,new_name,args.run)
    cmd(args.job, m_tag, matching, args.excluded, args.directory,rn_type,new_name,args.run)

if __name__ == "__main__":
    main()
