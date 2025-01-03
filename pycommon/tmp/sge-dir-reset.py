#!/gpfs/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import *

def run_commands(coms):
    for com in coms:
        os.system(com)
    return

def reset(package_name, dsave, save_files, Lrun):

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
        com = 'rm * '
        print(com)
        commands.append(com)
        com = 'rm -r *ampdb'
        print(com)
        commands.append(com)
        com = 'cp %s/* .' % dsave
        print(com)
        commands.append(com)

        if Lrun or yes_or_no("Did you save and will you run?"):
            run_commands(commands)
            
    return 0


def main():
    parser = argparse.ArgumentParser(description='reset directory with saving mode')
    parser.add_argument('package', choices=['amp'], help='choose packagename')
    parser.add_argument('-s', '--saving_dir', default='ini', help='save files to ./ini')
    parser.add_argument('-f', '--saving_files', nargs='*', help='select files to be saved')
    parser.add_argument('-r', '--run', action='store_true', help='run command if not just show')
    args = parser.parse_args()

    reset(args.package, args.saving_dir, args.saving_files, args.run)
    return 0

if __name__ == '__main__':
    main()
