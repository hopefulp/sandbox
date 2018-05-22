#!/home/jackjack5/epd/bin/python

import argparse
import subprocess as subp
import re
import os
import sys
from common import *

def cmd(job, match, lpattern, exceptions, dir_destination, new_name, run):

    if job == 'rm':
        print "use dir_clean.py"
        exit(1)
        
    pwd = os.getcwd()

    file_target=[]
    my_regex_l=[]
    n=0
    if match == 'p':
        f_list=get_files_prefix(lpattern, pwd)
    elif match == 's':
        f_list=get_files_suffix(lpattern, pwd)
    else:
        f_list=get_files_match(lpattern, pwd)

    if exceptions:
        for fex in exceptions:
            file_target.remove(fex)
    # run job for the file list
    #print("\n".join(file_target))
    if job == "ls":
        print("\n".join(f_list))
    elif job == "mv" :
        if not os.path.isdir(dir_destination):
            ret = subp.call(['mkdir', dir_destination])
            if ret == 1:
                sys.exit(1)
        for file in file_target:
            command = job + " " + file + " " + dir_destination
            if run:
                subp.Popen(command, shell=True, stdout=subp.PIPE)
            else:
                print command
    elif job == 'rename' or job == 'rn':
        if match == 'p':
            if len(prefixes) != 1:
                print "choose 1 prefix and new_name"
                exit(11)
                            
        else:
            n_pre=len(prefixes[0])
            for file in file_target:
                f_body = file[n_pre:]
                f_name  = prefixes[0] + f_body
                fn_name = new_name + f_body
                command  = "mv   %-30s %-30s" % (f_name, fn_name)
                if run:
                    subp.Popen(command, shell=True, stdout=subp.PIPE)
                else:
                    print command
    elif job == 'rm':
        for file in file_target:
            command = job + " " + file
            if run :
                subp.Popen(command, shell=True, stdout=subp.PIPE)
            else:
                print command
    
    if(job == 'rm' or job == 'mv' or job == 'rename') and run == False:
        print "use '-r' for run"
    return                    

#def run_bool(v):
#    if v.lower() in ('yes', '1', 'y', 't', 'true'):
#        return True

def main():
    parser = argparse.ArgumentParser(description='directory management')
    parser.add_argument( 'job', choices=['ls', 'mv', 'rm', 'rename', 'rn'],  help='shell command')
    #group=parser.add_mutually_exclusive_group()
    #parser.add_argument( '-p', '--prefix', nargs='*', help='prefix of filename')
    #parser.add_argument( '-m', '--match', nargs='*', help='find matching string')
    parser.add_argument('-m', '--matching', default='p', choices=['p', 's', 'm'], help='matching string with p-prefix, s-suffix, m-matching')
    parser.add_argument( '-s', '--string', nargs='*', help='list several searching word')
    
    parser.add_argument( '-e', '--exception', nargs='*', help='filename to be excluded')
    parser.add_argument( '-d', '--directory', type=str, default='tmppy', help='target directory to move files')
    parser.add_argument( '-rn', '--rename', help='rename filename')
    parser.add_argument( '-r', '--run', action='store_true', help='run or not-False')
    args = parser.parse_args()
    cmd(args.job, args.matching, args.string, args.exception, args.directory,args.rename,args.run)

if __name__ == "__main__":
    main()
