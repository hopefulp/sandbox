#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
import sys
import json
import subprocess
from subprocess import Popen, PIPE, STDOUT
from common     import get_dirs_prefix, yes_or_no, list2dict
from mod_incar  import modify_incar
from mod_poscar import fixedMD_POSCAR, pos2dirname, get_poscar
from mod_vas    import get_hostname, jg_poscar, jg_kpoints, jg_incar, jg_potcar, jg_link
from vas_qsub   import qsub_command

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

### O: Use this to comtain all the jobs
###                 1     2      3       4     5      6       7            8        9      10    11    12     13   14
def vasp_jobs(job, dirs, fixatom, kopt, iopt, ikw_opt, incar_kws, incar_list, Lrun, newdir, issue, np, xpart, nnode):
    pwd = os.getcwd()

    for odir in dirs:
        if not newdir:
            ndir = odir + job
        else:
            ndir = newdir
        com=[]
        ### 0: make a new dir
        #if os.path.isdir(ndir):
        #    print(f"{ndir} for {odir} exists: exits")
        #    sys.exit(1)
        #else:
        #    os.system(f'mkdir {ndir}')
        try: 
            subprocess.call(['mkdir', f'{ndir}'])     # str f'mkdir {ndir}' is not wokring
            print(f'{ndir} was made')
        except:
            print(f"can't make {ndir}")
            sys.exit(1)
        ### 1: POSCAR
        if job in jg_poscar:
            ### zpe: modify POSCAR
            if job == 'zpe':
                fixedMD_POSCAR(f"{odir}/CONTCAR", fixatom)
                print(f"{odir}/CONTCAR was modified to POSCAR")
                poscar = f'{pwd}/POSCAR'
            ### job == 'ini'
            else: 
                poscar = f'{odir}/POSCAR'
        ### other job uses CONTCAR
        else:
            poscar = f'{odir}/CONTCAR'
        os.system(f'cp {poscar} {ndir}/POSCAR')
        print(f"{poscar} was copied to {ndir}/POSCAR")

        ### 2: POTCAR
        if not job in jg_potcar:
            potcar = f'{odir}/POTCAR'
        else:
            # under construction
            pass
        os.system(f'cp {potcar} {ndir}')
        print(f"{potcar} was copied to {ndir}")

        ### 3: KPOINTS
        if not job in jg_kpoints:
            kpoints = f'{odir}/KPOINTS'
        else:
        #if vgroup == 1 or vgroup == 2:
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
        #elif vgroup == 3:
            elif job == 'zpe' or job == 'wav':
                kpoints = odir + '/KPOINTS'
            elif os.path.isfile(f'KPOINTS.{job}'):
                kpoints = 'KPOINTS.'+job
            elif os.path.isfile('KPOINTS'):
                kpoints = 'KPOINTS'
            else:
                kpoints = odir + "/KPOINTS"
        os.system(f'cp {kpoints} {ndir}/KPOINTS')
        print(f"{kpoints} was copied to {ndir}/KPOINTS")

        ### 4: INCAR
        #if (not ikw_opt or ikw_opt== 'm') and os.path.isfile(f"{odir}/INCAR"):
        #    incar = modify_incar(f"{odir}/INCAR", job)
        if not job in jg_incar:
            incar = f"{odir}/INCAR"
        else:
            if os.path.isfile(f"INCAR.{job}"):
                incar = f"INCAR.{job}"
            else:
                incar = make_incar(iopt, odir, job, ikw_opt, incar_kws, incar_list)
        os.system(f"cp {incar} {ndir}/INCAR")
        print(f"{incar} was copied to {ndir}/INCAR")

        ### 5: Link: 5: Link: 5: Link: 5: Link: 5: Link: make link for CHGCAR, WAVECAR
        if job in jg_link:
            os.chdir(ndir)
            if os.path.isfile(f'../{odir}/CHGCAR'):
                s = f"ln -s ../{odir}/CHGCAR ."
                os.system(s)
                print(f"CHGCAR is linked to {ndir}")
            else:
                print(f"there is no {ochgcar} for {job}")
                sys.exit(5) # exit number is file index 5 for CHGCAR, 6 for WAVECAR
            if job == 'pchg':
                if os.path.isfile(f'../{odir}/WAVECAR'):
                    s = f"ln -s ../{odir}/WAVECAR ."
                    os.system(s)
                    print(f"WAVECAR is linked to {ndir}")
                else:
                    print(f"there is no {owavecar} for {job}")
                    sys.exit(6)
            os.chdir(pwd)

        ### 6: Job submit: qsub
        qsub_command(ndir, X=xpart, nnode=nnode, np=np, issue=issue)
    return 0

def main():
    parser = argparse.ArgumentParser(description='remove files except initial files')
    parser.add_argument('-j', '--job', choices=['sp','incar',"dos","band","pchg","chg","chgw","md","cont","ini","zpe","mol","wav",'vdw','noD','opt','copt','mag','kisti'], help='inquire for each file ')
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
    parser.add_argument('-err', '--error', choices=['opt','mem','sim'], help="vasp error: converge, memory issue, sim for not to change INCAR")
    qsub = parser.add_argument_group(title='qsub')
    qsub.add_argument('-x', '--partition', type=int,  help='partition number in qsub')
    qsub.add_argument('-N', '--nnode',     type=int,  help='number of nodes in qsub')
    qsub.add_argument('-n', '--nproc',      help='nprocess in qsub')
    args = parser.parse_args()


    ### incar option
    if not args.ioption and 'opt' in args.job:
        ikw_option = 'ac'
    else:
        ikw_option = args.ioption

    ### copy initial job: POSCAR or CONTCAR
    ### cont + ok to change KPOINTS
    ### obtain job directories
    pwd = os.getcwd()
    if args.dirs:
        dirs = args.dirs
    elif args.prefix:
        dirs = get_dirs_prefix(pwd, args.prefix, excludes=args.exclude)
    else:
        print("Usage:: input old job dirs: -d ")
        sys.exit(1)

    vasp_jobs(args.job, dirs, args.fixed_atom, args.kopt, args.incar, ikw_option, args.incar_kws, args.incar_list, args.run,args.newdir,args.error, args.nproc, args.partition, args.nnode)

    return 0

if __name__ == '__main__':
    main()
