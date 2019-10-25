#!/usr/bin/python

import argparse
import os

def search_dir(pwd, dirs, i_files, f_more, L_rm):
    os.chdir(pwd)
    if f_more:
        i_files.extend(f_more)
        print i_files
    for dir in dirs:
        if dir != '.':
            os.chdir(dir)
            print("####.... enter %s directory" % dir)
        p_files = os.listdir('.')
        print p_files
        for file in p_files:
            if os.path.isdir(file):
                pass
            elif os.path.isfile(file):
                if file in i_files:
                    pass
                else:
                    cmd = 'rm  ' + file 
                    print cmd
                    if L_rm:
                        os.system(cmd)
        if dir != '.':                        
            os.chdir('..')                
            print "####.... exit %s directory" % dir

    if not L_rm:
        print "use -r to remove"
    return 0

def main():

    parser = argparse.ArgumentParser(description='remove files except initial files')
    parser.add_argument('-d', '--dir', nargs='*', default='.', help='select directorys')
    parser.add_argument('-f', '--files', nargs='*', help='files to be remode')
    parser.add_argument('-r', '--rm', action='store_true', help='remove or not option')
    args = parser.parse_args()

    ini_files = ['INCAR', 'POSCAR', 'POTCAR', 'KPOINTS']

    print args

    pwd = os.getcwd()
    os.chdir(pwd)
    dirs=[]
    if args.dir == '.':
        dirs.append(args.dir)
    elif args.dir == 'a':
        dirs = os.listdir(pwd)
    else:
        dirs = args.dir[:]
    print dirs

    search_dir(pwd, dirs, ini_files, args.files, args.rm)


    return 0

if __name__ == '__main__':
    main()
