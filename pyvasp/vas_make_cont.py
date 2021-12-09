#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
import json
from common import get_dirs_prefix, yes_or_no, list2dict
from mod_incar import modify_incar
from mod_poscar    import fixedMD_POSCAR, pos2dirname, get_poscar
import sys
from myvasp import get_hostname
from vas_qsub import qsub_command

### 1. for general job process
def vasp_jobs( job, dirs, prefix, exclude, fixatom, optin,opt_kp, Lrun, np, newdir):
    #print(f"{exclude}")
    pwd = os.getcwd()
    if prefix:
        dirs = get_dirs_prefix(pwd, prefix, excludes=exclude)

    for odir in dirs:
        if not newdir:
            ndir = odir + job
        else:
            ndir = newdir
        com=[]
        if os.path.isdir(ndir):
            print(f"{ndir} for {odir} exists")
            continue
        ### 0: make a new dir
        com.append(f'mkdir {ndir}')
        ### 1: POSCAR
        ### zpe: modify POSCAR
        if job == 'zpe':
            fixedMD_POSCAR(f"{odir}/CONTCAR", fixatom)
            print(f"{odir}/CONTCAR was modified to POSCAR")
            poscar = 'POSCAR'
        ### else: just copy CONTCAR
        else:
            poscar = odir + '/CONTCAR'
            print(f"{odir}/CONTCAR will be copied")
        com.append(f'cp {poscar} {ndir}/POSCAR')
        ### 2: POTCAR
        potcar = odir + '/POTCAR'
        com.append(f'cp {potcar} {ndir}/POTCAR')
        print(f"{potcar} was copied")
        ### 3: KPOINTS
        if job == 'zpe' or job == 'wav':
            kpoints = odir + '/KPOINTS'
        else:
            if os.path.isfile(f'KPOINTS.{job}'):
                kpoints = 'KPOINTS.'+job
            elif os.path.isfile('KPOINTS'):
                kpoints = 'KPOINTS'
        com.append(f'cp {kpoints} {ndir}/KPOINTS')
        print(f"{kpoints} was used")
        ### 4: INCAR
        if (not optin or optin== 'm') and os.path.isfile(f"{odir}/INCAR"):
            incar = modify_incar(f"{odir}/INCAR", job)
        elif optin == 'u' and os.path.isfile('INCAR.'+job):
            incar = 'INCAR.' + job
        else:
            incar = 'INCAR'
        com.append(f"cp {incar} {ndir}/INCAR")
        print(f"{incar} was used")

        ### make directory and copy
        for st in com:
            os.system(st)

        ### 5: More Extra files
        if job == 'band' or job == 'dos':
            os.chdir(ndir)
            s = f"ln -s ../{odir}/CHGCAR ."
            os.system(s)
            print(f"CHGCAR is linked to {ndir}")
            os.chdir(pwd)
        ### qsub depends on server
        s = qsub_command(ndir, np)
        print(s)
        if Lrun or yes_or_no("Will you run"):
            os.system(s)
    return 0


### 2 only incar is changed for jobs: vdw, 
def vasp_job_incar( job, dirs, prefix, exclude, fixatom, inc_option, Lrun, np, newdir, incar_kws, incar_list):
    '''
    in case only incar is changed
    job     vdw
    copy POTCAR KPOINTS CONTCAR change INCAR
        chg
    '''
    pwd = os.getcwd()
    if prefix:
        dirs = get_dirs_prefix(pwd, prefix, excludes=exclude)

    for odir in dirs:
        if newdir:
            ndir = newdir
        else:
            ndir = odir + job
        com=[]
        if os.path.isdir(ndir):
            print(f"{ndir} for {odir} exists")
            continue
        ### 0: make a new dir
        os.system(f'mkdir {ndir}')
        copyfiles = ['CONTCAR', 'KPOINTS', 'POTCAR']
        ### COPY 1,2,3 input files
        for f in copyfiles:
            if f == 'CONTCAR':
                if os.path.exists(f"{odir}/{f}") and os.stat(f"{odir}/{f}").st_size != 0:
                    com = f"cp {odir}/{f} {ndir}/POSCAR"
                else:
                    com = f"cp {odir}/POSCAR {ndir}"
            else:
                com = f"cp {odir}/{f} {ndir}/{f}"
            print(com)
            os.system(com)
        ### 4: INCAR
        if job == 'mag':
            dic = {'MAGMOM': inc_option[0]}
            opt='ac' # activate and change
        else:
            dic = {}
            opt = None
        if incar_kws:
            print(f"{incar_kws}")
            #kv = json.load(incar_kws)
            kws = list2dict(incar_kws)
            dic.update(kws)
        elif incar_list:
            dic = incar_list
            ### in case LDA: change POTCAR
            print(f"incar list {dic} and incar option {inc_option}")
            if 'GGA' in dic and 'o' in inc_option:
                s = "genpotcar.py -pp lda"
                os.chdir(ndir)
                os.system(s)
                os.chdir(pwd)
                print("POTCAR was changed with LDA")
            else:
                print("this is not working")
        if inc_option:
            opt = inc_option
        incar = modify_incar(f"{odir}/INCAR", job, dic=dic, opt=opt)
        print(f"{incar} was modified from {odir}/INCAR")
        com = f"cp {incar}  {ndir}/INCAR"
        print(f"{incar} was used")
        os.system(com)


    return 0

### 3:: only POSCAR or KPOINTS is changed for job ini & cont
def vasp_job_ini(job, dirs, poscar, newdir, Loptkp, Lrun):
    '''
    in case POSCAR or KPOINTS changes
    '''
    pwd = os.getcwd()

    odir = dirs[0]
    if newdir:
        ndir = newdir
    elif not poscar:
        if job == 'ini':
            ndir = odir+'n'
        elif job == 'cont':
            ndir = odir+'c'
        if Loptkp:
            ndir += 'kp'
    else:
        ndir = pos2dirname(poscar)
    if os.path.isdir(ndir):
        print(f"{ndir} for {odir} exists")
        sys.exit(1)
    ### 0: make a new dir
    os.system(f'mkdir {ndir}')
    copyfiles = ['INCAR', 'KPOINTS', 'POTCAR']
    ### COPY 1,2,3 input files
    for f in copyfiles:
        if f == 'KPOINTS' and Loptkp:
            print("prepare and use KPOINTS in wdir")
            com = f"cp KPOINTS {ndir}"
        else:
            com = f"cp {odir}/{f} {ndir}"
        print(f"{com}")
        os.system(com)
    ### 4: POSCAR
    if poscar:
        Poscar = poscar
    elif job == 'ini':
        Poscar = f"{odir}/POSCAR"
    elif job == 'cont':
        Poscar = f"{odir}/CONTCAR"
    com = f"cp {Poscar} {ndir}/POSCAR"
    print(f"{com}")
    os.system(com)
    ### only works for KISTI
    com = qsub_command(ndir)
    print(com)
    if Lrun:
        os.system(com)
    elif not poscar:
        print(f"made a new ini {ndir} directory")
    else:
        q = "will you run?"
        if yes_or_no(q):
            os.system(com)
    return 0


def main():
    parser = argparse.ArgumentParser(description='remove files except initial files')
    parser.add_argument('-j', '--job', choices=['sp','incar',"dos","band","pchg","chg","md","cont","ini","zpe","mol","wav",'vdw','noD','opt','copt','mag','kisti'], help='inquire for each file ')
    parser.add_argument('-d', '--dirs', nargs='+', help='select directories')
    parser.add_argument('-nd', '--newdir', help='select directories')
    parser.add_argument('-p', '--prefix', help='select directories using prefix')
    parser.add_argument('-ex', '--exclude', nargs='*', help='exclude if already exist')
    parser.add_argument('-a', '--fixed_atom', default='H', help='atom symbol to be fixed')
    parser.add_argument('-io', '--ioption', nargs='*', help='params a: append, c: change, o: out, r: reverse, u:use INCAR.job')
    #parser.add_argument('-id', '--incar_dict', type=json.loads, help='input dict from command line')
    parser.add_argument('-ikw', '--incar_kws', nargs='*', help='input key-value pairs in the list from command line')
    parser.add_argument('-il', '--incar_list', nargs='*', help='input list for comment out')
    parser.add_argument('-ok', '--optkpoints', action='store_true', help='make KPOINTS or copy KPOINTS.job')
    parser.add_argument('-s', '--poscar', help='incar POSCAR.name for job==ini')
    parser.add_argument('-r', '--run', action='store_true', help='Run without asking')
    parser.add_argument('-n', '--nproc', default=16, help='nprocess in qsub')
    args = parser.parse_args()


    ### incar option
    if not args.ioption and 'opt' in args.job:
        inc_option = 'ac'
    else:
        inc_option = args.ioption
    ### only INCAR changes in no-vdw -> vdw, sp -> opt, opt->sp
    incar_jobs = ['vdw','noD', 'opt','copt','mag', 'kisti','incar','sp','chg']
    ### copy initial job: POSCAR or CONTCAR
    ### cont + ok to change KPOINTS
    ini_jobs = ['ini', 'cont']
    ### others
    # band: INCAR + KPOINTS
    # dos:  INCAR + KPOINTS

    if args.job in incar_jobs:
        vasp_job_incar(args.job, args.dirs, args.prefix, args.exclude, args.fixed_atom, inc_option, args.run, args.nproc, args.newdir,args.incar_kws, args.incar_list)
    elif args.job in ini_jobs:
        vasp_job_ini( args.job, args.dirs, args.poscar, args.newdir, args.optkpoints, args.run)
    else:
        vasp_jobs(args.job, args.dirs, args.prefix, args.exclude, args.fixed_atom, inc_option,args.optkpoints, args.run, args.nproc, args.newdir )
    return 0

if __name__ == '__main__':
    main()
