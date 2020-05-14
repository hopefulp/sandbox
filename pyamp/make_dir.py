#!/home/joonho/anaconda3/bin/python

import argparse
import os
import sys

files = { 'md': ['OUTCAR', 'amp.amp']}

def jobs(d, work, job):
    pwd = os.getcwd()
    #directory = pwd + '/'+ d
    if os.path.isdir(d):
        print(f'there exist {d}')
        sys.exit(1)
    else:
        os.system(f'mkdir {d}')
        os.chdir(f'{d}')
    if work == 'amp':
        if job == 'md':
            for f in files['md']:
                forig = '../'+f
                if os.path.isfile(forig):
                    os.system(f'ln -s {forig} {f}')
                else:
                    print(f'there does not exist {forig}')
                    sys.exit(2)
                
    return 0

def main():
    parser = argparse.ArgumentParser(description="make directory with some auxiliary files  ")
    parser.add_argument( 'directory', default='tmp', help='make directory') 
    parser.add_argument('-w','--work', default='amp', choices=['amp'], help="work ")
    parser.add_argument('-j','--job', default='md', choices=['md'], help="md ")
    args = parser.parse_args()

    jobs(args.directory, args.work, args.job)

if __name__ == "__main__":
    main()
