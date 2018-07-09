#!/usr/bin/python

import argparse
import os

def vasp_jobs(dirs):
	host = os.getenv('HOST')

    if not dirs:
        print " run outsicde of vasp job directory"
        exit(0)
    elif dirs[0] == 'a':
        lists = os.listdir('.')
        for dir in lists:
            if not os.path.isdir(dir):
                continue
            os.chdir(dir)
            print 'now on '+dir
            if os.path.isfile('INCAR') and os.path.isfile('POSCAR') and not os.path.isfile('OUTCAR'):
                os.system('csh /qcfs/jackjack5/vasp/vaspenv.sh')
            os.chdir('..')
    else:
        for dir in dirs:
            if host == "" or host not in PPNs.keys():


def main():
	parser = argparse.ArgumentParser(description='execution of all vaspallallall')
    parser.add_argument('-d', '--dir', nargs='*', help='input directories, if not, present dir')
	args = parser.parse_args()
	vasp_jobs(args.dir)
	return



if __name__=='__main__':
	main()	


