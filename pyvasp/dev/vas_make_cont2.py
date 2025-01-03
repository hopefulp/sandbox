#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
import sys
import json
from common import get_dirs_prefix, yes_or_no, list2dict
from mod_incar import modify_incar
from mod_poscar    import fixedMD_POSCAR, pos2dirname, get_poscar
import sys
from envvasp import get_hostname
from vas_qsub import run_vasp

def change_incar(odir, ndir, job, incar_opt, incar_kws, incar_list):
    ### if incar_opt='u..', just use it
    if incar_opt and re.match('u', incar_opt):
        incar = 'INCAR.' + job
        print(f"{incar} is used")
    else:
        incaro = f"{odir}/INCAR"
        ### This should be < in case INCAR.job needs to be modified
        if job == 'mag':
            dic = {'MAGMOM': inc_option[0]}
            opt='ac' # activate and change
        else:
            dic = {}
            opt = None
        ### INCAR kw and list is exclusive option
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
        if incar_opt:
            opt = incar_opt
        incar = modify_incar(incaro, job, dic=dic, opt=opt)
        print(f"{incaro} was modified and {incar} was used")
    return incar

def make_incar(iopt, odir, job, ikw_opt, incar_kws, incar_list):
    print(f"{iopt}")
    if iopt:
        if iopt == 'job':
            job_incar = 'INCAR.' + job
            if os.path.isfile(job_incar):
                incar = job_incar
            else:
                print(f"There is no {incar} file to return")
                sys.exit(10)
        else:
            print(f"Note: iopt was defined but not 'job'")
    else:
        incar = change_incar(odir, job, ikw_opt, incar_kws, incar_list)
    return incar

### 1. INCAR and more files need to be modified
###                 1     2      3       4     5      6       7            8        9      10    11    12     13   14
def vasp_jobs_more(job, dirs, fixatom, kopt, iopt, ikw_opt, incar_kws, incar_list, Lrun, newdir, np, xpart, nnode, hmem):
    pwd = os.getcwd()

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
        if kopt:
            ### if kpoints file
            kfname = kopt
            kfsuff = f'KPOINTS.{kopt}'
            if os.path.isfile(kfname):
                kpoints = kfname
            ### if kpoints suffix
            elif os.path.isfile(kfsuff):
                kpoints = kfsuff
        ### use -j job
        else:
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
        #if (not ikw_opt or ikw_opt== 'm') and os.path.isfile(f"{odir}/INCAR"):
        #    incar = modify_incar(f"{odir}/INCAR", job)
        incar = make_incar(iopt, odir, job, ikw_opt, incar_kws, incar_list)
        st = f"cp {incar} {ndir}/INCAR"
        com.append(st)

        ### make directory and copy
        for st in com:
            print(f"{st}")
            os.system(st)

        ### 5: More Extra files
        if job == 'band' or job == 'dos':
            os.chdir(ndir)
            s = f"ln -s ../{odir}/CHGCAR ."
            os.system(s)
            print(f"CHGCAR is linked to {ndir}")
            os.chdir(pwd)
        ### qsub depends on server
        run_vasp(ndir, xpart, nnode, np, hmem)
    return 0


### 2 only incar is changed for jobs: vdw, 
def vasp_job_incar( job, dirs, fixatom, iopt, ikw_opt, Lrun,newdir, incar_kws, incar_list,np,xpart,nnode,hmem):
    '''
    in case only incar is changed
    job     vdw
    copy POTCAR KPOINTS CONTCAR change INCAR
        chg
    '''
    pwd = os.getcwd()

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
        incar = make_incar( iopt, odir, job, ikw_opt, incar_kws, incar_list)
        com = f"cp {incar} {ndir}/INCAR"
        print(f"{com}")
        os.system(com)

        ### qsub depends on server
        run_vasp(ndir, xpart, nnode, np, hmem)

    return 0

### 3:: only POSCAR or KPOINTS is changed for job ini & cont
def vasp_job_ini(job, dirs, poscar, newdir, Loptkp, Lrun, np, xpart, nnode, hmem):
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
    ### qsub depends on server
    run_vasp(ndir, xpart, nnode, np, hmem)

    return 0


def main():
    parser = argparse.ArgumentParser(description='remove files except initial files')
    parser.add_argument('-j', '--job', choices=['sp','incar',"dos","band","pchg","chg","md","cont","ini","zpe","mol","wav",'vdw','noD','opt','copt','mag','kisti'], help='inquire for each file ')
    dgroup = parser.add_mutually_exclusive_group()
    dgroup.add_argument('-d','-do', '--dirs', nargs='+', help='specify directories')
    dgroup.add_argument('-p', '--prefix', help='select directories using prefix')

    parser.add_argument('-ex', '--exclude', nargs='*', help='specify excluded dirs if already exist')
    parser.add_argument('-dn', '-nd', '--newdir', help='specify new dirname in case one job')

    parser.add_argument('-a', '--fixed_atom', default='H', help='atom symbol to be fixed')
    parser.add_argument('-i', '--incar', default='job', help='incar option:')
    parser.add_argument('-io', '--ioption', help='in the order: u:use INCAR.job,a:append,c:change,o:out,r:reverse')
    parser.add_argument('-ikw', '--incar_kws', nargs='*', help='input key-value pairs in the list from command line')
    parser.add_argument('-il', '--incar_list', nargs='*', help='input list for comment out')
    #kgroup = parser.add_mutually_exclusive_group('input kpoint option')
    parser.add_argument('-k', '--kopt', help='k option: fname, extension name, 3 values')
    #parser.add_argument('-k', '--optkpoints', action='store_true', help='make KPOINTS or copy KPOINTS.job')
    parser.add_argument('-s', '--poscar', help='incar POSCAR.name for job==ini')
    parser.add_argument('-r', '--run', action='store_true', help='Run without asking')
    qsub = parser.add_argument_group(title='qsub')
    qsub.add_argument('-x', '--partition', type=int,  help='partition number in qsub')
    qsub.add_argument('-N', '--nnode',     type=int,  help='number of nodes in qsub')
    qsub.add_argument('-n', '--nproc',      help='nprocess in qsub')
    qsub.add_argument('-m', '--hmem', action='store_true', help='in case large supercell, use half of memory')
    args = parser.parse_args()


    ### incar option
    if not args.ioption and 'opt' in args.job:
        ikw_option = 'ac'
    else:
        ikw_option = args.ioption
    ### only INCAR changes in no-vdw -> vdw, sp -> opt, opt->sp
    incar_jobs = ['vdw','noD', 'opt','copt','mag', 'kisti','incar','sp','chg']
    ### copy initial job: POSCAR or CONTCAR
    ### cont + ok to change KPOINTS
    ini_jobs = ['ini', 'cont']
    ### others
    # band: INCAR + KPOINTS
    # dos:  INCAR + KPOINTS
    ### obtain job directories
    pwd = os.getcwd()
    if args.dirs:
        dirs = args.dirs
    elif args.prefix:
        dirs = get_dirs_prefix(pwd, args.prefix, excludes=args.exclude)
    else:
        print("Usage:: input old job dirs: -d ")
        sys.exit(1)

    if args.job in incar_jobs:
        vasp_job_incar(args.job, dirs, args.fixed_atom, args.incar, ikw_option, args.run, args.newdir,args.incar_kws, args.incar_list,args.nproc,  args.partition, args.nnode, args.hmem)
    elif args.job in ini_jobs:
        vasp_job_ini( args.job, dirs, args.poscar, args.newdir, args.kopt, args.run, args.nproc,  args.partition, args.nnode, args.hmem)
    else:
        vasp_jobs_more(args.job, dirs, args.fixed_atom, args.kopt, args.incar, ikw_option, args.incar_kws, args.incar_list, args.run,args.newdir, args.nproc, args.partition, args.nnode, args.hmem )
    return 0

if __name__ == '__main__':
    main()
