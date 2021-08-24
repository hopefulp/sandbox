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


def make_vas_dir(old_dir, new_dir, job, incar, poscar, kpoints, potcar, Lvdw, Lquestion):

    pwd = os.getcwd()

    files=[]
    ### Define POSCAR, KPOINTS, INCAR
    if job == 'ini':
        poscar = 'POSCAR'
        kpoints = 'KPOINTS'
        incar = 'INCAR'
    else:
        poscar = 'CONTCAR'
        kpoints = 'IBZKPT'
        incar = 'INCAR'
        if job == 'zpe' or job == 'md' or job == 'chg':
            incar = 'INCAR.' + job 
            
        if os.path.isfile(incar):
            print(f"incar == {incar}, modify {incar}")
        else:
            print(f"No incar file {incar}")
            sys.exit(0)
    ### Define POTCAR
    if potcar:
        print("generate POTCAR")
    else:
        potcar = "POTCAR"
    ### Extra files for Continous Job
    #if job=="hybrid" or job=="md" or job=="dos":
    if job=="hybrid" or job=="md":
        files.append("WAVECAR")
    ### if vdW-DF
    if Lvdw:
        f_vdw = old_dir+'/'+"vdw_kernel.bindat"
        if os.path.isfile(f_vdw):
            files.append("vdw_kernel.bindat")

    if os.path.isdir(new_dir):
        print(f"overwrite {new_dir}")
    else:
        s = f"mkdir {new_dir}"
        os.system(s)
        print(f"directory {new_dir} was made")
    
    if Lquestion:
        ### Use the CONTCAR for POSCAR
        q = 'will you use CONTCAR?'
        if yes_or_no(q):
            pass
        elif yes_or_no('will you use POSCAR?'):
            poscar='POSCAR'    
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
        fo = old_dir +'/'+ poscar
        if os.path.isfile(fo):
            s = f"cp {fo} {new_dir}/POSCAR"
            os.system(s)
            print(f"{fo} will be copied to {new_dir}/POSCAR")
        else:
            print(f"There is not {fo}")

        ### POTCAR
        fo = old_dir + '/'  + potcar
        if os.path.isfile(fo):
            s = f"cp {fo} {new_dir}/POTCAR"
            os.system(s)
            print(f"{fo} will be copied to {new_dir}/POTCAR")
        else:
            print(f"There is not {fo}")
        ### KPOINTS
        fo = old_dir + '/'  + kpoints
        if job == 'zpe':
            fo = 'KPOINTS.gam'
        if os.path.isfile(fo):
            s = f"cp {fo} {new_dir}/KPOINTS"
            os.system(s)
            print(f"{fo} will be copied to {new_dir}/KPOINTS")
        else:
            print("make KPOINTS file ")
        ### INCAR :: copy INCAR or incar.key
        ### INCAR is normally different job so read on the work directory
        #fo = old_dir + '/'  + incar
        if os.path.isfile(incar) :
            s = f"cp {incar} {new_dir}/INCAR"
            os.system(s)
            print(f"{incar} will be  copied to {new_dir}/INCAR")
        else:
            print(f"There is not {incar}")
            
    if files:
        for f in files:
            fo = old_dir + '/'  + f
            s = f"cp {fo} {new_dir}"
            print(f"{fo} is copied to {new_dir}")
            os.system(s)
                
def main():
    global ini_dvasp, pwd
    parser = argparse.ArgumentParser(description='make directory for continous job from CONTCAR POTCAR KPOINTS and new INCAR')
    parser.add_argument('odir', help='copy from old dir')
    parser.add_argument('ndir', help='mkdir and cp')
    parser.add_argument('job', choices=["hybrid","dos","band","pchg","chg","md","cont","ini","zpe","mol"], help='inquire for each file')
    parser.add_argument('-s', '--poscar', help='use CONTCAR for 2nd job')
    parser.add_argument('-p', '--potcar', help='use the same POTCAR')
    parser.add_argument('-i', '--incar', default='INCAR', choices=["INCAR","incar.key"], help='first run make_incar.py then use incar.key')
    parser.add_argument('-k', '--kpoints', help='use the same KPOINTS for default')
    parser.add_argument('-q', '--question', action='store_true', help='inquire for each file')
    parser.add_argument('-v', '--vdw', action='store_true', help="in case vdW-DF,copy vdw_kernel.bindat")
    #parser.add_argument('-l', '--ksub', default='monk', choices=['monk','gamma','dos','band'], help='diverse k-point sampling')
    
    args = parser.parse_args()

    make_vas_dir(args.odir, args.ndir, args.job, args.incar, args.poscar, args.potcar, args.kpoints, args.vdw, args.question)
if __name__ == '__main__':
    main()
