#!/home/joonho/anaconda3/bin/python

import argparse
import os
import sys

amp_db = ['amp-fingerprint-primes.ampdb', 'amp-neighborlists.ampdb', 'amp-fingerprints.ampdb', 'OUTCAR']

files = { 'tr': amp_db, 'md': ['OUTCAR', 'amp.amp'], 'vasp': ['OUTCAR']}

def jobs(ndir, work, job, odir, Lcopy):
    pwd = os.getcwd()
    #directory = pwd + '/'+ ndir
    if os.path.isdir(ndir):
        print(f'there exist {ndir}')
        sys.exit(1)
    else:
        os.system(f'mkdir {ndir}')
        os.chdir(f'{ndir}')
    if work == 'amp':
        if job == 'md' or job == 'tr':
            for f in files[job]:
                #forig = '../'+f
                forig = pwd + '/' + f
                if os.path.isfile(forig):
                    if Lcopy:
                        os.system(f'cp {forig} {f}')
                    else:
                        os.system(f'ln -s {forig} {f}')
                else:
                    if f == 'amp.amp':
                        f = "amp-untrained-parameters.amp"
                        forig = pwd + '/' + f
                        if Lcopy:
                            os.system(f'cp {forig} {f}')
                        else:
                            os.system(f'cp {forig} {f}')
                    else:
                        print(f'there does not exist {forig}')
                        print(f"here is {os.getcwd()}")
                        sys.exit(11)
        elif job == 'vasp':
            if odir == None:
                print("input vasp directory")
                sys.exit(30)
            old_dir=pwd + '/' + odir
            for f in files[job]:
                forig = old_dir+'/'+f
                if os.path.isfile(forig):
                    os.system(f'cp {forig} .')
                else:
                    print(f'there does not exist {forig}')
                    sys.exit(31)
            
    return 0

def main():
    parser = argparse.ArgumentParser(description="make directory with some auxiliary files  ")
    parser.add_argument( 'directory', default='tmp', help='make directory') 
    parser.add_argument('-w','--work', default='amp', choices=['amp'], help="work ")
    parser.add_argument('-j','--job', default='md', choices=['tr','md','vasp'], help="md ")
    parser.add_argument('-od','--old_dir', help=" in case 'vasp', input old directory")
    parser.add_argument('-c', '--copy', action='store_true', help='make copy rather than ln')
    args = parser.parse_args()

    jobs(args.directory, args.work, args.job, args.old_dir, args.copy)

if __name__ == "__main__":
    main()
