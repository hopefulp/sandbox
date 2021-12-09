#!/home/joonho/anaconda3/bin/python
### versin 1.1 by J. Park
### 2018.4.2 makes input files by option -s(POSCAR) -p(POTCAR) -k(KPOINTS) -i(INCAR)
### incar is not ready
### 2019.10.25 update

import argparse
import os
import shutil
import re
from  envvasp import *
from common import *
from mod_poscar import get_poscar

def make_vas_d2d(odir, ndir, job, poscar, kpoints, potcar, incar):

    pwd = os.getcwd()
    if potcar and potcar=='lda':
        pp = 'lda'
    else:
        pp = 'gga'
    files=[]

    ### Extra files for Continous Job: CHGCAR, WAVECAR

    if os.path.isdir(ndir):
        print(f"overwrite {ndir}")
    else:
        s = f"mkdir {ndir}"
        os.system(s)
        print(f"directory {ndir} was made")
   
    ### 1. POSCAR
    get_poscar(poscar, job='new', sub=1)    #1 for CONTCAR
    s = f"cp POSCAR {ndir}"
    os.system(s)
    print(f"POSCAR  from {'wdir':>15} to {ndir}")
    ### 2. KPOINTS    
    if not kpoints:
        s = f"cp {odir:>15}/KPOINTS {ndir}"
    os.system(s)
    print(f"KPOINTS from {odir:>15} to {ndir}")
    ### 3. POTCAR
    if not potcar:
        s = f"cp {odir}/POTCAR {ndir}"
        os.system(s)
        print(f"POTCAR  from {odir:>15} to {ndir}")
    else:
        s = f"genpotcar.py -pp {pp}"
        os.chdir(ndir)
        os.system(s)
        os.chdir(pwd)
        print(f"POTCAR was generated in {ndir} with {pp}")
    ### INCAR :: copy INCAR or incar.key
    ### INCAR is normally different job so read on the work directory
    if incar:
        print("Use vas_make_cont.py")
        sys.exit(1)
    else:
        s = f"cp {odir}/INCAR {ndir}"
        os.system(s)
        print(f"INCAR   from {odir:>15} to {ndir}")
    return 0        
            
                
def main():
    global ini_dvasp, pwd
    parser = argparse.ArgumentParser(description='make directory from other directory ')
    parser.add_argument('odir', help='copy from old dir')
    parser.add_argument('ndir', help='mkdir and cp')
    parser.add_argument('-j', '--job', choices=['lda',"hybrid","md","mol"], help='inquire for each file')
    parser.add_argument('-s', '--poscar', help='copy odir/CONTCAR or input')
    parser.add_argument('-p', '--potcar', help='copy odir/potcaruse or make POTCAR')
    parser.add_argument('-i', '--incar',  help='use the same INCAR in d2d')
    parser.add_argument('-k', '--kpoints', help='copy odir/KPOINTS or make')
    
    args = parser.parse_args()
    make_vas_d2d(args.odir, args.ndir, args.job, args.poscar, args.kpoints, args.potcar, args.incar)
if __name__ == '__main__':
    main()
