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

def reset(package_name,command, dsave, save_files, Lrun):

    commands=[]
    if package_name == 'amp':
        '''
        if not os.path.isdir(dsave):
            print(" there is not dir %s" % dsave)
            os.mkdir(dsave)
            print("%s was created. Run again" % dsave)
            return

        if not os.listdir("./%s" % dsave):
            print("first save files to %s manually" % dsave)
            return
        '''
        if command == 'reset':
            com = 'rm * '
            print(com)
            commands.append(com)
            com = 'rm -r *ampdb'
            print(com)
            commands.append(com)
            com = 'cp %s/* .' % dsave
            print(com)
            commands.append(com)
        elif command == 'mkdir':
            
        if Lrun or yes_or_no("Did you save and will you run?"):
            run_commands(commands)
            
    return 0


def main():
    parser = argparse.ArgumentParser(description='reset directory with saving mode')
    parser.add_argument('-j', '--job', choices=['amp'], help='choose packagename')
    parser.add_argument('-c', '--command', choices=['reset','mkdir'], help='command=="mkdir": make directory and copy basic files\
    \n\t\t=="reset": save files and rerun at the same directory')
    parser.add_argument('-s', '--saving_dir', help='save files to ./dir')
    parser.add_argument('-f', '--saving_files', nargs='*', help='select files to be saved')
    parser.add_argument('-r', '--run', action='store_true', help='run command if not just show')
    args = parser.parse_args()

    if not args.job:
        print("input -j job such as 'amp'")
        sys.exit(0)

    reset(args.job, args.command, args.saving_dir, args.saving_files, args.run)
    return 0

if __name__ == '__main__':
    main()
