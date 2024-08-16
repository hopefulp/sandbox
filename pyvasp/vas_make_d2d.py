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
from mod_vas  import *

### vasp input order
### inputs = [args.poscar, args.kpoints, args.potcar, args.incar]
def make_vas_d2d(odir, ndir, job, inputs, files, np,xpart,nnode,hmem):

    vin_order = ['POSCAR', 'KPOINTS', 'POTCAR', 'INCAR']
    pwd = os.getcwd()
    if inputs[3] and inputs[3] =='lda':
        pp = 'lda'
    else:
        pp = 'gga'
  
    ### 1. POSCAR
    #get_poscar(poscar, job='new', sub=1)    #1 for CONTCAR

    vinputs = []
    opos = f"{odir}/POSCAR"
    
    for i, inp in enumerate(inputs):
        vin_old = f"{odir}/{vin_order[i]}"
        ### Use input
        if inp:
            vin = inp
        elif i == 0 and job == 'cont':    # for POSCAR
            if job == 'cont':
                vin = f"{odir}/CONTCAR"
                if os.path.isfile(vin):
                    print(f"{vin} for {vin_order[i]}")
        ### Use old directory
        elif os.path.isfile(vin_old):
                vin = vin_old 
                print(f"{vin} for {vin_order[i]}")
        else:
            vin = None
            print(f"there is no {vin_order[i]}")
            
        vinputs.append(vin)

    ### Make ndir
    if os.path.isdir(ndir):
        print(f"overwrite {ndir}")
    else:
        s = f"mkdir {ndir}"
        os.system(s)
        print(f"directory {ndir} was made")

    ### copy files
    for i, vin in enumerate(vinputs):
        if vin:
            s = f"cp {vin} {ndir}/{vin_order[i]}"
            os.system(s)
            print(f"{s}")
        else:
            if i == 2:
                s = f"genpotcar.py -pp {pp}"
                os.chdir(ndir)
                os.system(s)
                os.chdir(pwd)
                print(f"POTCAR was generated in {ndir} with {pp}")
    ### Extra files for Continous Job: CHGCAR, WAVECAR
    if files:
        for f in files:
            s = f"cp -P {odir}/{f} {ndir}"
            os.system(s)
            print(f"{odir}/{f} was copied to {ndir}")
    ### copy subdirectiry with only POSCAR
    if job == 'neb':
        subfiles = os.listdir(odir)
        #print(f"old dir: {subfiles}")
        for f in subfiles:
            if os.path.isdir(f"{odir}/{f}"):
                os.system(f"mkdir {ndir}/{f}")
                os.system(f"cp {odir}/{f}/POSCAR {ndir}/{f}")
            elif re.search('POSCAR', f):
                os.system(f"cp {odir}/{f} {ndir}")
    ### run?
    run_vasp(ndir, xpart, nnode, np, hmem)

    return 0        
            
                
def main():
    global ini_dvasp, pwd
    parser = argparse.ArgumentParser(description='make directory from other directory ')
    parser.add_argument('odir', help='copy from old dir')
    parser.add_argument('ndir', help='mkdir and cp')
    parser.add_argument('-j', '--job', choices=['lda','hybrid','md','mol','kp','ini','neb','cont'], help='inquire for each file')
    parser.add_argument('-s', '--poscar', help='copy POSCAR')
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
    inputs = [args.poscar, args.kpoints, args.potcar, args.incar]
    make_vas_d2d(args.odir, args.ndir, args.job, inputs, args.files, args.partition, args.nnode, args.nproc, args.hmem)
if __name__ == '__main__':
    main()
