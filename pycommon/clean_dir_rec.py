#!/usr/bin/python

import argparse
import os

def search_dir(job, dir, foption, lfile_exc, exe_tag):
    os.chdir(dir)
    print("####.... enter %s directory" % dir)
    p_files = os.listdir('.')
    #print("%s" % type(foption))
    for file in p_files:
        if not os.path.isdir(file):
            if file in lfile_exc:
                pass
            else:
                if isinstance(foption, list):
                    if file in foption:
                        cmd = job + ' ' + file 
                        print cmd
                        if exe_tag:
                            os.system(cmd)
                            print("%s was removed" % file)
                    else:
                        pass
                elif isinstance(foption, float):
                    fsize = os.path.getsize(file)
                    fsize /= 1024.0 * 1024.0
                    #print("%f %f" % (foption, fsize))
                    if fsize > foption:
                        cmd = job + ' ' + file
                        print cmd
                        if exe_tag:
                            os.system(cmd)
                            print("%s was removed" % file)
                    else:
                        pass
        # if directory
        else:
            search_dir(job, file, foption, lfile_exc, exe_tag)
    os.chdir('..')                
    print "####....... exit %s directory" % dir
def main():

    parser = argparse.ArgumentParser(description='remove files recursively u. size or filename')
    parser.add_argument('job', choices=['rm'], help='job is only for remove')
    parser.add_argument('dir', nargs='*', help='select directorys')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-f', '--files', nargs='*', help='files to be remode')
    group.add_argument('-s', '--size', type=float, default=10.0, help='delete file of large size, default=10M')
    parser.add_argument('-x', '--xfiles', nargs='*', default=['INCAR', 'POSCAR', 'POTCAR', 'KPOINTS', 'IBZKPT', 'CONTCAR'], help='excluded files')
    parser.add_argument('-e', '--exe', action='store_true', help='remove or not option')
    args = parser.parse_args()

    job=args.job

    print args

    pwd = os.getcwd()
    os.chdir(pwd)

    del_files = ['WAVECAR', 'CHGCAR', 'CHG']

    if args.files:
        foption = args.files
    else:
        foption = args.size

    for dir in args.dir:
        search_dir(job, dir, foption, args.xfiles, args.exe)

    return 0

if __name__ == '__main__':
    main()
