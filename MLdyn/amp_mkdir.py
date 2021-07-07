#!/home/joonho/anaconda3/bin/python

import argparse
import os
import sys
import amp_util

amp_db = ['amp-fingerprint-primes.ampdb', 'amp-neighborlists.ampdb', 'amp-fingerprints.ampdb', 'OUTCAR']
### 'db' for what?
files = { 'tr': amp_db, 'md': ['OUTCAR'], 'vasp': ['OUTCAR'], 'db': amp_db, 'te': amp_db, 'des':['OUTCAR'] }

def jobs(ndir, work, job, odir, Lcopy):
    pwd = os.getcwd()
    #directory = pwd + '/'+ ndir
    if os.path.isdir(ndir):
        print(f'there exist {ndir}')
        sys.exit(1)
    #### make new directory 
    else:
        os.system(f'mkdir {ndir}')
        os.chdir(f'{ndir}')
    if work == 'amp':
        if job == 'md' or job == 'tr' or job == 'db' or job == 'te' or job == 'des':
            for f1 in files[job]:
                f = pwd+'/'+f1
                if os.path.isfile(f) or os.path.isdir(f):
                    if Lcopy:
                        os.system(f'cp ../{f1} .')
                    else:
                        os.system(f'ln -s ../{f1} .')   # link should be used from link file location
                else:
                    print(f"can't find {f1}")
                    sys.exit(0)
            if job == 'md' or job == 'te':
                f = amp_util.get_amppotname()
                if f:
                    if Lcopy:                    
                        os.system(f'cp {f} {ndir}')
                    else:
                        os.system(f'ln -s ../{f} {ndir}/{f}')
                else:
                    forig = pwd + '/amp-untrained-parameters.amp'
                    f = 'amp-untrained-parameters.amp'
                    os.system(f'ln -s {forig} {f}')
                    print(f"No amp file in {pwd}")
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
            
            
            
    os.chdir(f'{pwd}')
    return 0

def main():
    parser = argparse.ArgumentParser(description="make directory with some auxiliary files  ")
    parser.add_argument( 'directory', default='tmp', help='make directory') 
    parser.add_argument('-w','--work', default='amp', choices=['amp'], help="work ")
    parser.add_argument('-j','--job', default='md', choices=['des','tr','md','db','vasp','te'], help="des: different descriptor test ")
    parser.add_argument('-od','--old_dir', help=" in case 'vasp', input old directory")
    parser.add_argument('-c', '--copy', action='store_true', help='make copy rather than ln')
    args = parser.parse_args()

    jobs(args.directory, args.work, args.job, args.old_dir, args.copy)

if __name__ == "__main__":
    main()
