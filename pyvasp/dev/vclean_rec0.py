#!/usr/bin/python

import argparse
import os

def search_dir(job, dir, lfile_del, lfile_exc, exe_tag):
    os.chdir(dir)
    print(("####.... enter %s directory" % dir))
    p_files = os.listdir('.')
    for file in p_files:
        if not os.path.isdir(file):
            if file in lfile_exc:
                pass
            elif file in lfile_del:
                cmd = job + ' ' + file 
                print(cmd)
                if exe_tag:
                    os.system(cmd)
            else:
                pass
        # if directory
        else:
            search_dir(job, file, lfile_del, lfile_exc, exe_tag)
    os.chdir('..')                
    print("####.... exit %s directory" % dir)
def main():

    parser = argparse.ArgumentParser(description='remove files recursively in subdirectories')
    parser.add_argument('job', choices=['rm'], help='job is only for remove')
    parser.add_argument('dir', nargs='*', help='select directorys')
    parser.add_argument('-f', '--files', nargs='+', default=['WAVECAR','CHGCAR','CHG'], help='files to be remode')
    parser.add_argument('-x', '--xfiles', nargs='*', default=['INCAR', 'POSCAR', 'POTCAR', 'KPOINTS', 'IBZKPT', 'CONTCAR'], help='excluded files')
    parser.add_argument('-e', '--exe', action='store_true', help='remove or not option')
    args = parser.parse_args()

    job=args.job

    print(args)

    pwd = os.getcwd()
    os.chdir(pwd)

    for dir in args.dir:
        search_dir(job, dir, args.files, args.xfiles, args.exe)

    return 0

if __name__ == '__main__':
    main()
