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
    parser.add_argument('odir', help='copy from old dir')
    parser.add_argument('ndir', help='mkdir and cp')
    parser.add_argument('-q', '--question', action='store_true', help='inquire for each file')
    parser.add_argument('-j', '--job', choices=["hybrid","dos","band","pchg","md","cont","ini"], help='inquire for each file')
    parser.add_argument('-v', '--vdw', action='store_true', help="in case vdW-DF,copy vdw_kernel.bindat")
    parser.add_argument('-s', '--poscar', default="CONTCAR", help='use CONTCAR for 2nd job')
    parser.add_argument('-p', '--potcar', default="POTCAR" , help='use the same POTCAR')
    parser.add_argument('-k', '--kpoints', default="IBZKPT", help='use the same KPOINTS for default')
    #parser.add_argument('-l', '--ksub', default='monk', choices=['monk','gamma','dos','band'], help='diverse k-point sampling')
    parser.add_argument('-i', '--incar', default='INCAR', choices=["INCAR","incar.key"], help='first run make_incar.py then use incar.key')
    args = parser.parse_args()

    pwd = os.getcwd()
    old_dir = args.odir
    new_dir = args.ndir 

    files=[]
    if args.job:
        if args.job=="hybrid" or args.job=="md" or args.job=="dos":
            files.append("WAVECAR")
    f_vdw = old_dir+'/'+"vdw_kernel.bindat"
    if args.vdw or os.path.isfile(f_vdw):
        files.append("vdw_kernel.bindat")

    if os.path.isdir(new_dir):
        print(f"overwrite {new_dir}")
    else:
        s = f"mkdir {new_dir}"
        os.system(s)
        print(f"directory {new_dir} was made")
    
    if args.question:
        ### Use the CONTCAR for POSCAR
        q = 'will you use CONTCAR?'
        if yes_or_no(q):
            pass
        elif yes_or_no('will you use POSCAR?'):
            args.poscar='POSCAR'    
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
        fo = old_dir +'/'+ args.poscar
        if not os.path.isfile(fo):
            print(f"There is not {fo}")
            fo = old_dir + '/' + "POSCAR"
        s = f"cp {fo} {new_dir}/POSCAR"
        os.system(s)
        print(f"{fo} was copied to {new_dir}/POSCAR")

        ### POTCAR
        fo = old_dir + '/'  + args.potcar
        if os.path.isfile(fo):
            s = f"cp {fo} {new_dir}"
            os.system(s)
            print(f"{fo} was copied to {new_dir}")
        else:
            print(f"There is not {fo}")
            sys.exit(3)
        ### KPOINTS
        fo = old_dir + '/'  + args.kpoints
        if os.path.isfile(fo):
            s = f"cp {fo} {new_dir}/KPOINTS"
            os.system(s)
            print(f"{fo} was copied to {new_dir}/KPOINTS")
        else:
            fo = old_dir + '/'  + "KPOINTS"
            if os.path.isfile(fo):
                s = f"cp {fo} {new_dir}/KPOINTS"
                os.system(s)
                print(f"{fo} was copied to {new_dir}/KPOINTS")
            else:
                print("make KPOINTS ")
                sys.exit(4)
        ### INCAR :: copy INCAR or incar.key       
        fo = old_dir + '/'  + args.incar
        if os.path.isfile(fo) :
            s = f"cp {fo} {new_dir}"
            os.system(s)
            print(f"{fo} was copied to {new_dir}")
        else:
            print(f"There is not {fo}")
        print("Modify INCAR for different Job")
            
    if files:
        for f in files:
            fo = old_dir + '/'  + f
            s = f"cp {fo} {new_dir}"
            print(f"{fo} is copied to {new_dir}")
            os.system(s)
                
if __name__ == '__main__':
    main()
