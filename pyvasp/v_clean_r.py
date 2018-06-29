#!/home/jackjack5/epd/bin/python

import argparse
import os

def search_dir(job, dir, lfile_save, lfile_del):
    os.chdir(dir)
    print("####.... enter %s directory" % dir)
    p_files = os.listdir('.')
    for file in p_files:
        if not os.path.isdir(file):
            if file in lfile_save:
                pass
            elif file in lfile_del:
                cmd = job + ' ' + file 
                print cmd
                os.system(cmd)
            else:
                pass
        # if directory
        else:
            search_dir(job, file, lfile_save, lfile_del)
    os.chdir('..')                
    print "####.... exit %s directory" % dir
def main():

    parser = argparse.ArgumentParser(description='remove files except initial files')
    parser.add_argument('job', choices=['rm'], help='job is only for remove')
    parser.add_argument('dir', nargs='*', help='select directorys')
    parser.add_argument('-f', '--files', nargs='+', default=['WAVECAR','CHGCAR','CHG'], help='files to be remode')
    args = parser.parse_args()

    job=args.job

    ini_files = ['INCAR', 'POSCAR', 'POTCAR', 'KPOINTS', 'IBZKPT', 'CONTCAR']

    print args

    pwd = os.getcwd()
    os.chdir(pwd)

    for dir in args.dir:
        search_dir(job, dir, ini_files, args.files)


    return 0

if __name__ == '__main__':
    main()
