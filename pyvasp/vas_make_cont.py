#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import get_dirs_prefix, yes_or_no
from mod_incar import modify_incar
from mod_poscar    import fixedMD_POSCAR, pos2dirname, get_poscar
import sys
from myvasp import get_hostname
from vas_qsub import qsub_command

incar_job = ['vdw', 'opt', 'mag']

def vasp_jobs( job, dirs, prefix, exclude, fixatom, opt_incar, Lrun, np, newdir):
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
        if job == 'zpe':
            modify_incar(f"{odir}/INCAR", job)
            incar = 'INCAR'
            print(f"{incar} was modified from {odir}/INCAR")
        else:
            if (not opt_incar or opt_incar== 'm') and os.path.isfile(f"{odir}/INCAR"):
                incar = modify_incar(f"{odir}/INCAR", job)
            elif opt_incar == 'c' and os.path.isfile('INCAR.'+job):
                incar = 'INCAR.' + job
            else:
                incar = 'INCAR'
        com.append(f"cp {incar} {ndir}/INCAR")
        print(f"{incar} was used")

        ### make directory and copy
        for st in com:
            os.system(st)

        ### 5: More Extra files
        if job == 'band':
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



def vasp_job_incar( job, dirs, prefix, exclude, fixatom, option, Lrun, np, newdir):
    '''
    in case only incar is changed
    job     vdw
    '''
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
        os.system(f'mkdir {ndir}')
        copyfiles = ['CONTCAR', 'KPOINTS', 'POTCAR']
        ### COPY 1,2,3 input files
        for f in copyfiles:
            if f == 'CONTCAR':
                com = f"cp {odir}/{f} {ndir}/POSCAR"
            else:
                com = f"cp {odir}/{f} {ndir}/{f}"
            os.system(com)
        ### 4: INCAR
        if job in incar_job:
            if job == 'mag':
                dic = {'MAGMOM': option[0]}
                opt='ac' # activate and change
            else:
                dic = None
                opt = None
            incar = modify_incar(f"{odir}/INCAR", job, dic=dic, opt=opt)
            print(f"{incar} was modified from {odir}/INCAR")
        com = f"cp {incar}  {ndir}/INCAR"
        print(f"{incar} was used")
        os.system(com)
    return 0


def vasp_job_pos(dirs, poscar, newdir, Lrun):
    '''
    in case POSCAR or KPOINTS changes
    '''
    pwd = os.getcwd()

    odir = dirs[0]
    if newdir:
        ndir = newdir
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
        com = f"cp {odir}/{f} {ndir}"
        os.system(com)
    ### 4: POSCAR
    com = f"cp {poscar} {ndir}/POSCAR"
    os.system(com)
    ### only works for KISTI
    com = qsub_command(ndir)
    print(com)
    if Lrun:
        os.system(com)
    else:
        q = "will you run?"
        if yes_or_no(q):
            os.system(com)
    return 0


def main():
    parser = argparse.ArgumentParser(description='remove files except initial files')
    parser.add_argument('-j', '--job', choices=["dos","band","pchg","chg","md","cont","ini","zpe","mol","wav",'vdw','opt','mag'], help='inquire for each file')
    parser.add_argument('-d', '--dirs', nargs='+', help='select directories')
    parser.add_argument('-nd', '--newdir', help='select directories')
    parser.add_argument('-p', '--prefix', help='select directories using prefix')
    parser.add_argument('-ex', '--exclude', nargs='*', help='exclude if already exist')
    parser.add_argument('-a', '--fixed_atom', default='H', help='atom symbol to be fixed')
    parser.add_argument('-o', '--option', nargs='*', help='incar option, file option for ini ')
    parser.add_argument('-s', '--poscar', help='incar POSCAR.name for job==ini')
    parser.add_argument('-r', '--run', action='store_true', help='Run without asking')
    parser.add_argument('-n', '--nproc', default=16, help='nprocess in qsub')
    args = parser.parse_args()

    ### only INCAR changes in no-vdw -> vdw, sp -> opt,
    if args.job in incar_job :
        vasp_job_incar(args.job, args.dirs, args.prefix, args.exclude, args.fixed_atom, args.option, args.run, args.nproc, args.newdir )
    elif args.job == 'ini':
        if args.poscar:
            vasp_job_pos( args.dirs, args.poscar, args.newdir, args.run)
        else:
            print("input -s POSCAR.name")
            sys.exit(1)
    else:
        vasp_jobs(args.job, args.dirs, args.prefix, args.exclude, args.fixed_atom, args.option, args.run, args.nproc, args.newdir )
    return 0

if __name__ == '__main__':
    main()
