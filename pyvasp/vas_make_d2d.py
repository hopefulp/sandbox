#!/home/joonho/anaconda3/bin/python
### versin 1.1 by J. Park
### 2018.4.2 makes input files by option -s(POSCAR) -p(POTCAR) -k(KPOINTS) -i(INCAR)
### incar is not ready
### 2019.10.25 update

import argparse
import os
import shutil
import re
from common import *
from mod_poscar import get_poscar
from vas_qsub import run_vasp
from vas_env  import *


def make_vas_d2d(odir, ndir, job, poscar, kpoints, potcar, incar, files, np,xpart,nnode,hmem):

    pwd = os.getcwd()
    if potcar and potcar=='lda':
        pp = 'lda'
    else:
        pp = 'gga'


    if os.path.isdir(ndir):
        print(f"overwrite {ndir}")
    else:
        s = f"mkdir {ndir}"
        os.system(s)
        print(f"directory {ndir} was made")
   
    ### 1. POSCAR
    #get_poscar(poscar, job='new', sub=1)    #1 for CONTCAR
    if poscar == 'p':
        pos = f"{odir}/POSCAR"
    elif poscar == 'c':
        pos = f"{odir}/CONTCAR"
    elif poscar == 'w':
        pos = f"{pwd}/POSCAR"
    else:
        if os.path.isfile(poscar):
            pos = poscar
        else:
            print("there is no poscar")
    s = f"cp {pos} {ndir}/POSCAR"
    os.system(s)
    print(f"POSCAR  from {pos:>15} to {ndir}")
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
    if incar:
        inc = incar
    else:
        inc = f"{odir}/INCAR"
    s = f"cp {inc} {ndir}/INCAR"
    os.system(s)
    print(f"INCAR   from {inc:>15} to {ndir}")
    
    ### Extra files for Continous Job: CHGCAR, WAVECAR
    if files:
        for f in files:
            s = f"cp -P {odir}/{f} {ndir}"
            os.system(s)
            print(f"{odir}/{f} was copied to {ndir}")
    ### run?
    run_vasp(ndir, xpart, nnode, np, hmem)

    return 0        
            
                
def main():
    global ini_dvasp, pwd
    parser = argparse.ArgumentParser(description='make directory from other directory ')
    parser.add_argument('odir', help='copy from old dir')
    parser.add_argument('ndir', help='mkdir and cp')
    parser.add_argument('-j', '--job', choices=['lda','hybrid','md','mol','kp'], help='inquire for each file')
    parser.add_argument('-s', '--poscar', default='p', help='p[POSCAR],c[CONTCAR],w[wdir/POSCAR],input poscar')
    parser.add_argument('-p', '--potcar', help='copy odir/potcaruse or make POTCAR')
    parser.add_argument('-i', '--incar',  help='use the same INCAR in d2d')
    parser.add_argument('-k', '--kpoints', help='copy odir/KPOINTS or make')
    parser.add_argument('-f', '--files', nargs='*', help='copy more files')
    qsub = parser.add_argument_group(title='qsub')
    qsub.add_argument('-x', '--partition',  help='partition number in qsub')
    qsub.add_argument('-N', '--nnode',      help='number of nodes in qsub')
    qsub.add_argument('-n', '--nproc',      help='nprocess in qsub')
    qsub.add_argument('-m', '--hmem', action='store_true', help='in case large supercell, use half of memory')
      
    args = parser.parse_args()
    make_vas_d2d(args.odir, args.ndir, args.job, args.poscar, args.kpoints, args.potcar, args.incar, args.files, args.partition, args.nnode, args.nproc, args.hmem)
if __name__ == '__main__':
    main()
