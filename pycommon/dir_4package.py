#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
import sys
from common import *

def run_commands(coms):
    for com in coms:
        os.system(com)
    return

def reset(package_name,command, dirname, save_files, Lrun):

    l_commands=[]
    if package_name == 'amp':
        '''
        if not os.path.isdir(dirname):
            print(" there is not dir %s" % dirname)
            os.mkdir(dirname)
            print("%s was created. Run again" % dirname)
            return

        if not os.listdir("./%s" % dirname):
            print("first save files to %s manually" % dirname)
            return
        '''
        if not command:
            print("use -c for command")
        elif command == 'reset':
            com = 'rm * '
            print(com)
            l_commands.append(com)
            com = 'rm -r *ampdb'
            print(com)
            l_commands.append(com)
            com = 'cp ini/* .'  
            print(com)
            l_commands.append(com)
        elif command == 'rerun':
            i=0
            dname = command + i
            #os.path.is
        elif command == 'mkdir':
            l_commands.append(f'mkdir {dirname}')
            l_commands.append(f'cp *ampdb *extxyz {dirname}')

            
        if l_commands != []:
            #print(f"l_commands {l_commands}")
            if Lrun or yes_or_no("Did you save and will you run?"):
                run_commands(l_commands)
            
    return 0


def main():
    parser = argparse.ArgumentParser(description='reset directory with saving mode')
    parser.add_argument('-j', '--job', choices=['amp'], help='choose packagename')
    parser.add_argument('-c', '--command', choices=['reset','mkdir','rerun'], help='command=="mkdir": make directory and copy basic files\
    \n\t\t=="reset": save files and rerun at the same directory')
    parser.add_argument('-d', '--dirname', help='directory for saving files or making news directory')
    parser.add_argument('-f', '--saving_files', nargs='*', help='select files to be saved')
    parser.add_argument('-r', '--run', action='store_true', help='run command if not just show')
    args = parser.parse_args()

    if not args.job:
        print("input -j job such as 'amp'")
        sys.exit(0)

    reset(args.job, args.command, args.dirname, args.saving_files, args.run)
    return 0

if __name__ == '__main__':
    main()
