#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
import sys
import json
import subprocess
from subprocess import Popen, PIPE, STDOUT
from common     import get_dirfiles, yes_or_no, list2dict
from libincar  import modify_incar_byjob, modify_incar_bykv
from libposcar import modify_POSCAR, pos2dirname, get_poscar
from libvas    import get_hostname, jg_poscar, jg_kpoints, jg_incar, jg_potcar, jg_link
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
            kws = list2dict(incar_kws)
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
def vasp_cont_1dir(job, odir, ndir, incar_kws, kopt, Lrun, option, np, xpart, nnode):
    '''
    job     default = 'cont'
            depending on job, some files will be changed
    odir
    ndir
    kopt
    incar_kws   dict with keys of 'i'|'a'|'d' for add or delete
                'i': if INCAR is defined
                'a': dict of k:w to be add or active
                'd': delete key list

    '''
    pwd = os.getcwd()

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
    if job in jg_poscar:    # use POSCAR
        ### zpe: modify POSCAR
        if job == 'zpe':
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

    ###### 4: INCAR
    #if (not ikw_opt or ikw_opt== 'm') and os.path.isfile(f"{odir}/INCAR"):
    #    incar = modify_incar_byjob(f"{odir}/INCAR", job)
    ### define INCAR to be used
    ### 4.1 Use it, if defined
    if 'i' in incar_kws.keys():
        if os.path.isfile(f"{incar_kws['i']}"):
            incar = f"{incar_kws['i']}"
        elif os.path.isdir(f"{incar_kws['i']}") and os.path.isfile(f"{incar_kws['i']}/INCAR"):
            incar = f"{incar_kws['i']}/INCAR"
    ### 4.2 try INCAR.job 
    elif os.path.isfile(f"INCAR.{job}"):
        incar = f"INCAR.{job}"
    ### 4.3 Use odir/INCAR for job = cont, ...
    elif not job in jg_incar:
        incar = f"{odir}/INCAR"
    ### if incar needs to be modified: kw for active, remove for comment out
    if 'a' in incar_kws or 'd' in incar_kws:
        ### make incar.new
        incar_o = incar
        if 'a' in incar_kws:
            incar = modify_incar_bykv(incar_o, incar_kws['a'], mode='m')   # return output filename
        if 'd' in incar_kws:
            incar = modify_incar_bykv(incar_o, incar_kws['d'], mode='e')
        print(f"{incar_o} was modified to {incar}")
    else:
        ### change incar by job

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
        if job == 'pchg' or job == 'cont':
            if os.path.isfile(f'../{odir}/WAVECAR'):
                s = f"ln -s ../{odir}/WAVECAR ."
                os.system(s)
                print(f"WAVECAR is linked to {ndir}")
            else:
                print(f"there is no {owavecar} for {job}")
                sys.exit(6)
        os.chdir(pwd)

    ### 6: Job submit: qsub
    qsub_command(ndir, X=xpart, nnode=nnode, np=np, option=option)
    return 0

def main():
    parser = argparse.ArgumentParser(description='How to make a continuous job dir')
    ### job
    parser.add_argument('-j', '--job', default='cont', choices=['sp','cont','incar','dos','band','pchg','chg','chgw','md','ini','kp','zpe','mol','wav','vdw','noD','opt','copt','mag','kisti'], help='inquire for each file ')
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
    goutput.add_argument('-a', '--fixed_atom', help='atom symbol to be fixed')      # default='H'
    ### modify 4 files
    ### INCAR
    parser.add_argument('-i', '--incar', help='specify incar file or dir: m for modify')
    parser.add_argument('-ia', '--incar_add', nargs='*', help='input key-value pairs in the list from command line')
    parser.add_argument('-id', '--incar_del', nargs='*', help='input list for comment out')
    #kgroup = parser.add_mutually_exclusive_group('input kpoint option')
    parser.add_argument('-k', '--kopt', help='k option: fname, extension name, 3 values')
    #parser.add_argument('-k', '--optkpoints', action='store_true', help='make KPOINTS or copy KPOINTS.job')
    parser.add_argument('-pos', '--poscar', help='incar POSCAR.name for job==ini')
    parser.add_argument('-r', '--run', action='store_true', help='Run without asking')
    parser.add_argument('-o', '--option', choices=['opt','mem','sim','long'], help="vasp error: converge, memory issue, sim for not to change INCAR")
    qsub = parser.add_argument_group(title='qsub')
    qsub.add_argument('-x', '--partition', type=int,            help='partition number in qsub')
    qsub.add_argument('-N', '--nnode',     type=int,            help='number of nodes in qsub')
    qsub.add_argument('-np', '--nproc',                          help='nprocess in qsub')
    parser.add_argument('-u', '--usage',   action='store_true', help='print usage')
    args = parser.parse_args()

    pwd = os.getcwd()
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
                ndir = odir+"fixed"+args.fixed_atom
            else:
                ndir = odir+f'{args.job}'
            new_dirs.append(ndir)
            #i += 1
    ### treat INCAR
    incar_kw={}
    if args.incar:
        incar_kw['i'] = args.incar          # value = string for input INCAR file
    if args.incar_add:
        incar_add = list2dict(args.incar_add)    # value = dict for k-w pair 
        incar_kw['a'] = incar_add
    if args.incar_del:
        incar_kw['d'] = args.incar_del      # value = list
        
    ### check directories
    print(f'{old_dirs} {new_dirs}')
    for odir, ndir in zip(old_dirs, new_dirs):
        print(f'run {odir} to {ndir}')
        ### call single jobs
        ###         1         2         3        4          5            6           7              8           9           10 
        vasp_cont_1dir(args.job, odir, ndir, incar_kw, args.kopt, args.run,  args.option, args.nproc, args.partition, args.nnode)

    return 0

if __name__ == '__main__':
    main()
