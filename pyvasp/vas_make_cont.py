#!/home/joonho/anaconda3/bin/python

import argparse
import os
from common import get_dirs_prefix, yes_or_no
from mod_vas    import fixedMD_POSCAR
import sys
from myvasp import get_hostname
from vas_qsub import qsub_command
from incar_mod import modify_incar

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
            if opt_incar:
                incar = modify_incar(f"{odir}/INCAR", job)
            elif os.path.isfile(f"INCAR.{job}"):
                incar = 'INCAR.' + job
            elif os.path.isfile('INCAR'):
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

def main():
    parser = argparse.ArgumentParser(description='remove files except initial files')
    parser.add_argument('-j', '--job', choices=["dos","band","pchg","chg","md","cont","ini","zpe","mol","wav"], help='inquire for each file')
    parser.add_argument('-d', '--dirs', nargs='+', help='select directories')
    parser.add_argument('-nd', '--newdir', help='select directories')
    parser.add_argument('-p', '--prefix', help='select directories using prefix')
    parser.add_argument('-ex', '--exclude', nargs='*', help='exclude if already exist')
    parser.add_argument('-a', '--fixed_atom', default='H', help='atom symbol to be fixed')
    parser.add_argument('-i', '--incar', help='Future incar option ')
    parser.add_argument('-r', '--run', action='store_true', help='Run without asking')
    parser.add_argument('-n', '--nproc', default=16, help='nprocess in qsub')
    args = parser.parse_args()

    vasp_jobs(args.job, args.dirs, args.prefix, args.exclude, args.fixed_atom, args.incar, args.run, args.nproc, args.newdir )
    return 0

if __name__ == '__main__':
    main()
