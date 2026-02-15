#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
import sys
import json
import subprocess
from subprocess import Popen, PIPE, STDOUT
from common     import get_dirfiles, yes_or_no 
from libstr     import li2dic
from libincar   import modify_incar_byjob, modify_incar_bykv, add_inckv_bysubjob
from libposcar  import modify_POSCAR, pos2dirname, get_poscar
from libvas     import get_hostname, jg_poscar, jg_kpoints, jg_incar, jg_potcar, jg_linkw, jg_linkc
from libkpoints import mod_kpoints
from vas_qsub   import qsub_command

def change_incar(odir, ndir, job, incar_opt, incar_kws, incar_remove):
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
            kws = li2dic(incar_kws)
            dic.update(kws)

        elif incar_remove:
            dic = incar_remove
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
        incar = modify_incar_byjob(incaro, job, dic=dic, opt=opt)
        print(f"{incaro} was modified and {incar} was used")
    return incar

def make_incar(iopt, odir, job, ikw_opt, incar_kws, incar_remove):
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
        incar = change_incar(odir, job, ikw_opt, incar_kws, incar_remove)
    return incar

     

### O: Use this to comtain all the jobs
###            1     2      3       4          5      6      7       8      9    10 
def vasp_cont_1dir(job, subjob, odir, ndir, incar, incopt, kopt, Lrun, option, np, xpart, nnode):
    '''
    job     default = 'cont'
            depending on job, some files will be changed
    odir
    ndir
    kopt
    incar       use input incar file
    incopt    if value == None, delete key
    '''
    pwd = os.getcwd()

    vjob = job
    if subjob:
        vjob += subjob
    ### treat INCAR
    if incopt:
        incopt = li2dic(incopt)    # value = dict for k-w pair 
    ### 0: make a new dir
    try: 
        subprocess.call(['mkdir', f'{ndir}'])     # str f'mkdir {ndir}' is not wokring
        print(f'{ndir} was made')
    except:
        print(f"can't make {ndir}")
        sys.exit(10)
    ### 0: subdir
    odir_full = f'{pwd}/{odir}'
    print(odir_full)
    ### 1: POSCAR
    if vjob in jg_poscar:    # use POSCAR
        ### zpe: modify POSCAR
        if vjob == 'zpe':
            modify_POSCAR(f"{odir}/CONTCAR", fixatom)
            print(f"{odir}/CONTCAR was modified to POSCAR")
            poscar = f'{pwd}/POSCAR'
        ### job == 'ini'
        else: 
            poscar = f'{odir}/POSCAR'
    ### other job uses CONTCAR
    else:
        ### subdirs is generated in NEB
        print("in case CONTCAR")
        if 'subdirs' in locals():
            for subdir in subdirs:
                poscar = f'{odir}/{subdir}/CONTCAR'
                os.system(f'cp {poscar} {ndir}/{subdir}/POSCAR')
        else:
            poscar = f'{odir}/CONTCAR'
    os.system(f'cp {poscar} {ndir}/POSCAR')
    print(f"{poscar} was copied to {ndir}/POSCAR in case no subdirectories")

    ### 2: POTCAR
    if not vjob in jg_potcar:
        potcar = f'{odir}/POTCAR'
    else:
        # under construction
        pass
    os.system(f'cp {potcar} {ndir}')
    print(f"{potcar} was copied to {ndir}")

    ### 3: KPOINTS
    ### copy old dir: zpe, wav, pchg
    ### change KPONTS in jg_kpoints: 'dos','kp'
    old_kpoints = odir + '/KPOINTS'
    
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
    elif os.path.isfile(f'KPOINTS.{vjob}'):
        kpoints = 'KPOINTS.'+vjob
    elif os.path.isfile('KPOINTS'):
        kpoints = 'KPOINTS'
    
    if 'kpoints' in locals():
        os.system(f'cp {kpoints} {ndir}/KPOINTS')
        print(f"{kpoints} was copied to {ndir}/KPOINTS")
    else:
        print(f"{old_kpoints} is modified to KPOINTS.new")
        if vjob in jg_kpoints:
            new_kpoints = 'KPOINTS.new'
            mod_kpoints(old_kpoints, outf=new_kpoints, kmulti=3)
            os.system(f'cp {new_kpoints} {ndir}/KPOINTS')
            print(f"{new_kpoints} was copied to {ndir}/KPOINTS")
        else:
            print(f"no rule for making KPOINTS")
            sys.exit(3)

    ###### 4: INCAR
    ### 4.1 Use it, if defined
    if incar: 
        if os.path.isfile(f"{incar}"):
            pass
        elif os.path.isdir(f"{incar}") and os.path.isfile(f"{incar}/INCAR"):
            incar = f"{incar}/INCAR"
    ### 4.2 try INCAR.job 
    elif os.path.isfile(f"INCAR.{vjob}"):
        incar = f"INCAR.{vjob}"
    #elif os.path.isfile(f"INCAR.{job}"):
    #    incar = f"INCAR.{job}"
    ### 4.3 Use odir/INCAR for job = cont, ...
    #elif not job in jg_incar:
    else:
        incar = f"{odir}/INCAR"
    ### if incar needs to be modified: kw for active, remove for comment out
    ### link job + subjob
    print(f"job {job} fullname-job {vjob} in {__file__}")
    if not os.path.isfile(f"INCAR.{vjob}"):
        if vjob in jg_incar:
            print(f"modify {incar} byjob")
            modify_incar_byjob(incar, vjob, outf=f'INCAR.{vjob}')
            incar = f'INCAR.{vjob}'
            #print(f"modify {incar} by job {vjob}")
        else:
            add_inckv_bysubjob(job, subjob, incopt)
    
    ### additional change 
    if incopt:
        ### make modify the present INCAR file
        print(f"modify {incar} by input -io")
        modify_incar_bykv(incar, incopt, outf='INCAR.new', mode='m')
        incar = 'INCAR.new'
        print(f"modify {incar} by option {incopt}")

    os.system(f"cp {incar} {ndir}/INCAR")
    print(f"{incar} was copied to {ndir}/INCAR")
    ### change into INCAR.newjob
    #if incar == 'INCAR.new':
    #    os.system(f"mv {incar} INCAR.{job}")
    #    print(f"{incar} was renamed to INCAR.{job}")

    ### 5: Link: 5: Link: 5: Link: 5: Link: 5: Link: make link for CHGCAR, WAVECAR
    if vjob in jg_linkw or vjob in jg_linkc:
        os.chdir(ndir)
        ### select CHGCAR or WAVECAR
        if vjob in jg_linkw and os.path.isfile(f'../{odir}/WAVECAR'):
            s = f"ln -s ../{odir}/WAVECAR ."
            os.system(s)
            print(f"WAVECAR is linked to {ndir}")
        else:
            print("can't link to WAVECAR")
        if vjob in jg_linkc and os.path.isfile(f'../{odir}/CHGCAR'):
            s = f"cp ../{odir}/CHGCAR ."
            os.system(s)
            print(f"CHGCAR is linked to {ndir}")
        else:
            print("can't link to CHGCAR")
        os.chdir(pwd)

    ### 6: Job submit: qsub
    qsub_command(ndir, X=xpart, nnode=nnode, np=np, option=option)
    return 0

def main():
    '''
    input   spw[a2]     after scf write WAVECAR, CHGCAR(2), PROCAR(a-default)
    '''
    parser = argparse.ArgumentParser(description='How to make a continuous job dir')
    ### job
    parser.add_argument('-j', '--job', default='cont', choices=['sp','cont','dos','band','chg','pchg','md','mdnve','nnff','nnffnve','ini','kp','zpe','mol','wav','vdw','noD','opt','copt','mag','kisti'], help='inquire for each file ')
    parser.add_argument('-sj', '--subjob', choices=['w','wlp','w2','sp', 'cool', 'heat','quench','B','kisti'],\
             help='sp for fake and others for md, for pchg, B for Bader')
        ### old directory selection
    gdirectory = parser.add_mutually_exclusive_group()
    gdirectory.add_argument('-d','-do', '--dirs', nargs='+', help='specify directories')
    gdirectory.add_argument('-p', '--prefix', help='select directories using prefix')
    ### exception for directory selection
    parser.add_argument('-ex', '--exclude', nargs='*', help='specify excluded dirs if already exist')
    ### job directory naming
    goutput = parser.add_mutually_exclusive_group()
    goutput.add_argument('-n', '--newdirs', nargs='+', help='specify new dirname in case one job')
    goutput.add_argument('-s', '--suffix', help='specify suffix of new directory')
    goutput.add_argument('-fa', '--fixed_atom', help='atom symbol to be fixed')      # default='H'
    ### modify 4 files
    ### INCAR
    parser.add_argument('-i', '--incar', help='specify incar file or dir: m for modify')
    parser.add_argument('-io', '--incar_option', nargs='*', help='input key-value pairs in the list from command line, if value=None, del')
    #kgroup = parser.add_mutually_exclusive_group('input kpoint option')
    parser.add_argument('-k', '--kopt', help='k option: fname, extension name, 3 values')
    #parser.add_argument('-k', '--optkpoints', action='store_true', help='make KPOINTS or copy KPOINTS.job')
    parser.add_argument('-pos', '--poscar', help='incar POSCAR.name for job==ini')
    parser.add_argument('-r', '--run', action='store_true', help='Run without asking')
    qsub = parser.add_argument_group(title='qsub')
    qsub.add_argument('-x', '--partition', type=int,            help='partition number in qsub')
    qsub.add_argument('-N', '--nnode',     type=int,            help='number of nodes in qsub')
    qsub.add_argument('-np', '--nproc',                          help='nprocess in qsub')
    qsub.add_argument('-o', '--option', choices=['opt','mem','sim','long','ml'], help="vasp error: converge, memory issue, sim for not to change INCAR")
    parser.add_argument('-u', '--usage',   action='store_true', help='print usage')
    args = parser.parse_args()

    pwd = os.getcwd()

    #print(f"after parse_args: fixed atoms {args.fixed_atom}") # error: -al conceived by -a l

    ### get old dirnames list
    if args.dirs:
        old_dirs = args.dirs
    elif args.prefix:
        old_dirs = get_dirs(pwd, prefix=args.prefix, excludes=args.exclude)
    else:
        print("Usage:: input old job dirs: -d ")
        sys.exit(1)

    ### get new dirnames
    new_dirs=[]
    #i=0
    if args.newdirs:
        new_dirs.extend(args.newdirs)
    else:
        for odir in old_dirs:
            if args.suffix:
                ndir = odir+args.suffix
            elif args.fixed_atom:
                print(f"fixed atoms {args.fixed_atom}")
                ndir = odir+"fixed"+args.fixed_atom
            else:
                ndir = odir+f'{args.job}'
                if args.subjob:
                    ndir += f'{args.subjob}'
            new_dirs.append(ndir)
            #i += 1
        
    ### check directories
    print(f'{old_dirs} {new_dirs}')
    for odir, ndir in zip(old_dirs, new_dirs):
        print(f'run {odir} to {ndir}')
        ### call single jobs
        ###         1         2         3        4          5            6           7              8           9           10 
        vasp_cont_1dir(args.job, args.subjob, odir, ndir, args.incar, args.incar_option, args.kopt, args.run,  args.option, args.nproc, args.partition, args.nnode)

    return 0

if __name__ == '__main__':
    main()
