#!/home/joonho/anaconda3/bin/python
### versin 1.1 by J. Park
### 2018.4.2 makes input files by option -s(POSCAR) -p(POTCAR) -k(KPOINTS) -i(INCAR)
### incar is not ready
### 2019.10.25 update

import argparse
import os
import shutil
import re
from  myvasp import *
from common import *

def main():
    global ini_dvasp, pwd
    parser = argparse.ArgumentParser(description='make directory for continous job from CONTCAR POTCAR KPOINTS and new INCAR')
    parser.add_argument('dir', help='mkdir and cp')
    parser.add_argument('-q', '--question', action='store_true', help='inquire for each file')
    parser.add_argument('-j', '--job', choices=["hybrid","dos","band","pchg"], help='inquire for each file')
    parser.add_argument('-s', '--poscar', default="CONTCAR", help='use CONTCAR for 2nd job')
    parser.add_argument('-p', '--potcar', default="POTCAR" , help='use the same POTCAR')
    parser.add_argument('-k', '--kpoints', default="IBZKPT", help='use the same KPOINTS for default')
    #parser.add_argument('-l', '--ksub', default='monk', choices=['monk','gamma','dos','band'], help='diverse k-point sampling')
    parser.add_argument('-i', '--incar', default='INCAR', choices=["INCAR","incar.key"], help='first run make_incar.py then use incar.key')
    parser.add_argument('-o', '--old_dir', help='copy from old dir')
    args = parser.parse_args()

    pwd = os.getcwd()

    files=[]
    if args.job:
        if args.job == "hybrid":
            files.append("WAVECAR")

    if not args.dir:
        print("input directory:")
        sys.exit(1)
    elif os.path.isdir(args.dir):
        print(f"overwrite {args.dir}")
    else:
        s = f"mkdir {args.dir}"
        os.system(s)
        print(f"directory {args.dir} was made")
    
    if args.old_dir:
        pre_dir =  args.old_dir+"/"
    else:
        pre_dir = "./"

    s = f"cp {pre_dir}*.sh {args.dir}"
    os.system(s)
        
    if args.question:
        ### Use the CONTCAR for POSCAR
        q = 'will you use CONTCAR?'
        if yes_or_no(q):
            pass
        else:
            print("What else?")
            sys.exit(0)
        ### Use the same POTCAR
        ### Use the same KPOINTS
        q = 'will you use the same KPOINTS?'
        if yes_or_no(q):
            pass
        else:
            q = 'input nummber of kpoints: [gamma|3 digits]'
            kp_in = get_answers(q)
            if re.match("g", kp_in, re.IGNORECASE) :
                method = "gamma"
                kps = "1  1  1"
            else:
                lkp = kp_in.split()
                if len(lkp) == 3:
                    method = 'MH'
                    kps = kp_in
                    print('default is MH')
                else:
                    print("input error for KPOINTS")
                    exit(11)
            print(kps, method)
            make_kpoints(kps, method)
    ### without -q            
    else:
        ### CONTCAR
        fo = pre_dir + args.poscar
        if os.path.isfile(fo):
            s = f"cp {fo} {args.dir}/POSCAR"
            os.system(s)
        else:
            print(f"There is not {fo}")
            sys.exit(2)
        ### POTCAR
        fo = pre_dir + args.potcar
        if os.path.isfile(fo):
            s = f"cp {fo} {args.dir}"
            os.system(s)
        else:
            print(f"There is not {fo}")
            sys.exit(3)
        ### KPOINTS
        fo = pre_dir + args.kpoints
        if os.path.isfile(fo):
            s = f"cp {fo} {args.dir}/KPOINTS"
            os.system(s)
        else:
            print("make KPOINTS ")
            sys.exit(4)
        ### INCAR :: copy INCAR or incar.key       
        fo = pre_dir + args.incar
        if os.path.isfile(fo) :
            s = f"cp {fo} {args.dir}"
            os.system(s)
        else:
            print(f"There is not {fo}")
            sys.exit(5)
    if files:
        for f in files:
            fo = pre_dir + f
            s = f"cp {fo} {args.dir}"
            os.system(s)
                
if __name__ == '__main__':
    main()
