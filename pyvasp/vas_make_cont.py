#!/home/joonho/anaconda3/bin/python

import argparse
import os
from common import get_dirs_prefix, yes_or_no
from mod_vas    import fixedMD_POSCAR


def vasp_jobs( job, dirs, prefix, exclude, fixatom, Lincar, Lrun, np):
    #print(f"{exclude}")
    pwd = os.getcwd()
    if prefix:
        dirs = get_dirs_prefix(pwd, prefix, excludes=exclude)

    if job == 'zpe':
        for odir in dirs:
            ndir = odir + 'zpe'
            commands=[]
            commands.append(f'mkdir {ndir}')
            commands.append(f'cp {odir}/POTCAR {ndir}/')
            commands.append(f'cp {odir}/KPOINTS {ndir}/')
            
            fixedMD_POSCAR(f"{odir}/CONTCAR", fixatom)
            commands.append(f'cp POSCAR {ndir}')
            
            ### make INCAR or use INCAR.zpe
            if Lincar:
                modify_INCAR(f"{odir}/INCAR", job)
                commands.append( f'cp INCAR {ndir}')
            else:
                commands.append(f"cp INCAR.{job} {ndir}/INCAR")
                print(f"INCAR.{job} is copied")
            
            ### qsub
            commands.append(f"qsub -N {ndir} -pe numa {np} -v np={np} -v dir={ndir} -v vas=gam $SB/pypbs/sge_vasp_exe.csh")
            for command in commands:
                print(command)
            if Lrun or yes_or_no("Will you run"):
                for command in commands:
                    os.system(command)
    return 0

def main():

    parser = argparse.ArgumentParser(description='remove files except initial files')
    parser.add_argument('-j', '--job', choices=["hybrid","dos","band","pchg","chg","md","cont","ini","zpe","mol"], help='inquire for each file')
    parser.add_argument('-d', '--dirs', nargs='*', help='select directories')
    parser.add_argument('-p', '--prefix', help='select directories using prefix')
    parser.add_argument('-ex', '--exclude', nargs='*', help='exclude if already exist')
    parser.add_argument('-a', '--fixed_atom', default='H', help='atom symbol to be fixed')
    parser.add_argument('-i', '--Lincar', action='store_true', help='T: make incar, F: copy INCAR.job')
    parser.add_argument('-r', '--run', action='store_true', help='Run without asking')
    parser.add_argument('-n', '--nproc', default=16, help='nprocess in qsub')
    args = parser.parse_args()

    vasp_jobs(args.job, args.dirs, args.prefix, args.exclude, args.fixed_atom, args.Lincar, args.run, args.nproc )
    return 0

if __name__ == '__main__':
    main()
