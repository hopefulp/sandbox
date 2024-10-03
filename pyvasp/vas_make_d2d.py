#!/home/joonho/anaconda3/bin/python
### versin 1.1 by J. Park
### 2018.4.2 makes input files by option -s(POSCAR) -p(POTCAR) -k(KPOINTS) -i(INCAR)
### incar is not ready
### 2019.10.25 update

import argparse
import os, sys
import shutil
import re
from common import *
from mod_poscar import get_poscar
from vas_qsub import qsub_command
from mod_vas  import *
from incar_change import change_incar_byjob

### vasp input order
### inputs = [args.poscar, args.kpoints, args.potcar, args.incar]
def make_vas_d2d(odir, ndir, job, inputs, files, qx, qN, qn, option=None, vasp_exe=None, lkisti=None, Lrun=None):

    if not os.path.isdir(odir):
        print(f"there is not {odir} directory, then, stop")
        sys.exit(11)

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
    if 'neb' in job:
        subfiles = os.listdir(odir)
        #print(f"old dir: {subfiles}")
        ### f dirname: 00 01 02 ...
        for f in subfiles:
            if os.path.isdir(f"{odir}/{f}"):
                os.system(f"mkdir {ndir}/{f}")
                if 'cont' in job:
                    pos = f"{odir}/{f}/CONTCAR"
                    if os.path.exists(pos):
                    ### Do not use WAVECAR, CHGCAR
                        #os.chdir(f"{ndir}/{f}")
                        ### WAVECAR is linked -> change INCAR
                        #os.system(f"ln -s ../../{odir}/{f}/WAVECAR .")
                        #os.chdir(f"{pwd}")
                        pass
                    else:
                        pos = f"{odir}/{f}/POSCAR"
                    os.system(f"cp {pos} {ndir}/{f}/POSCAR")
            ### copy ini and fin POSCAR in neb directory
            elif re.search('POSCAR', f):
                os.system(f"cp {odir}/{f} {ndir}")
        os.chdir(ndir)
        change_incar_byjob('INCAR', 'cont', outf='INCAR')
        os.chdir(pwd)
    ### run?

    #run_vasp(ndir, xpart, nnode, np, option)
    if get_hostname()=='pt' and (not qx or not qN):
        qx, qN = get_queue_pt(qx=qx, opt='long')
    s = qsub_command(ndir,X=qx,nnode=qN, np=qn, option=option, vasp_exe=vasp_exe, lkisti=lkisti, Lrun=Lrun)
    return 0        
            
                
def main():
    parser = argparse.ArgumentParser(description='make directory from other directory ')
    parser.add_argument('odir', help='copy from old dir')
    parser.add_argument('ndir', help='mkdir and cp')
    parser.add_argument('-j', '--job', choices=['lda','hybrid','md','mol','kp','ini','cont','neb','nebcont'], help='inquire for each file')
    parser.add_argument('-s', '--poscar', help='copy POSCAR')
    parser.add_argument('-p', '--potcar', help='copy odir/potcaruse or make POTCAR')
    parser.add_argument('-i', '--incar',  help='use the same INCAR in d2d')
    parser.add_argument('-k', '--kpoints', help='copy odir/KPOINTS or make')
    parser.add_argument('-f', '--files', nargs='*', help='copy more files')
    qsub = parser.add_argument_group(title='qsub')
    qsub.add_argument('-x', '--partition',  help='partition number in qsub')
    qsub.add_argument('-N', '--nnode',      help='number of nodes in qsub')
    qsub.add_argument('-n', '--nproc',      help='nprocess in qsub')
    #qsub.add_argument('-m', '--hmem', action='store_true', help='in case large supercell, use half of memory')
    qsub.add_argument('-o', '--option', choices=['long', 'mem'], help='option for qsub command line input')
    parser.add_argument('-u', '--usage',   action='store_true', help='print usage')

    args = parser.parse_args()
    inputs = [args.poscar, args.kpoints, args.potcar, args.incar]

    if args.usage:
        print(f"neb::\
                \n\t{__file__.split('/')[-1]} {args.odir} {args.ndir}  -j [neb|nebcont] ")
        sys.exit(1)

    #make_vas_d2d(args.odir, args.ndir, args.job, inputs, args.files, args.partition, args.nnode, args.nproc, args.hmem)
    make_vas_d2d(args.odir, args.ndir, args.job, inputs, args.files, args.partition, args.nnode, args.nproc, option=args.option)
if __name__ == '__main__':
    main()
