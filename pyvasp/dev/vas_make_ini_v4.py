#!/home/joonho/anaconda3/bin/python
### versin 1.1 by J. Park
### 2018.4.2 makes input files by option -s(POSCAR) -p(POTCAR) -k(KPOINTS) -i(INCAR)
### incar is not ready
### 2019.10.25 update
### 2021.05.07 update for EE

import argparse
import os
import shutil
import re
from common     import *
from vas_qsub   import get_queue_pt, qsub_command
from libvas    import *
from mod_poscar import get_poscar, pos2dirname

home = os.environ['HOME']
hostname = get_hostname()
pseudo_pot={'new':'Pot-new', 'potpaw-pbe-new':'Pot-new', 'old':'pot-old', 'potpaw-pbe-old':'pot-old'}
global pwd, ini_dvasp

#    copy {poscar} to POSCAR at cwd
def get_potcar(pot,atoms):
    potfiles=[]
    potdir = ini_dvasp + '/' + pseudo_pot[pot]
    files = os.listdir(potdir)
    #print 'atoms are ', atoms
    for atom in atoms:
        d_tag = 'No'
        for f1 in files:
            st_list = f1.split('_')     # atomic potcar file name is delimited by '_'
            if len(st_list) >= 1:
                if st_list[-1] == atom:
                    d_tag = 'Yes'
            if d_tag == 'Yes':                  
                potfiles.append(f1)
                break
        if d_tag == 'No':
            print("atom %s is not detected" % atom)
            exit(32)
    comm = 'cat '            
    for potfile in potfiles:
        fname = potdir + '/' + potfile
        comm += fname + ' '
    comm += ' > POTCAR'
    os.system(comm)
    print('POTCAR is combined from atoms ', atoms)

    return 0    

def get_incar(ifile):
    dic = {}
    make_incar(dic, ifile)
    print('INCAR is made from %s' % ifile)

    return 0

def make_vasp_dir(job, poscar, apotcar, hpp_list, kpoints, opt_incar, allprepared, dirname, iofile, atoms, issue, Lrun,Lmkdir,qx,qN,qn):
    global ini_dvasp, pwd
    ### 0. obtain default vasp repository
    ini_dvasp = get_vasp_repository()
    pwd = os.getcwd()
    files2copy=[]
    ### Now this is running
    ### 1. get POSCAR: make dirname using poscar
    if not poscar:
        q = 'will you make POSCAR? '
        if yes_or_no(q):
            q = 'input file: '
            poscar = get_answers(q)
    else:
        ### cp input poscar to 'POSCAR'
        #print(f"{__name__}::{poscar}")
        get_poscar(poscar)
        if not dirname:
            dirname = pos2dirname(poscar)
    ### use filename as it is
    files2copy.append('POSCAR')
    ### 2. get POTCAR will be made from scratch
    # POTCAR will be made after POSCAR moves to job directory
    ### 3. get KPOINTS
    q = 'will you make KPOINTS?'
    kpointjob = 'KPOINTS.' + job
    q2 = f'will you use {kpointjob}?'
    if kpoints:
        if len(kpoints) == 3: 
            method = "MH"
        elif len(kpoints) == 1 and kpoints[0] == 'g':
            kpoints = "1  1  1"
            method = "gamma"
        make_kpoints(kpoints, method)
    elif allprepared:
        print(f"KPOINTS in cwd will be copied to {dirname}")
    elif os.path.isfile(kpointjob) and yes_or_no(q2):
        os.system(f'cp {kpointjob} KPOINTS')
    elif yes_or_no(q):
        q = 'input nummber of kpoints: [gamma|3 digits such as "4 4 1" ]? '
        kp_in = get_answers(q)
        if re.match("g", kp_in, re.IGNORECASE) :
            method = "gamma"
            kps = "1  1  1"
        else:
            lkp = kp_in.split()
            if len(lkp) == 3:
                method = 'MH'
                kps = kp_in
                print('default is MH')
            else:
                print(f"input error for KPOINTS: {lkp} of length {len(lkp)}")
                exit(11)
        print(kps, method)
        make_kpoints(kps, method)
    else:
        print("Use wdir KPOINTS")
    files2copy.append('KPOINTS')
    ### 4. get INCAR :: use make_incar.py
    incar_repo=f"{ini_dvasp}/INCAR.{job}"
    incar_job=f"INCAR.{job}"
    if opt_incar:
        ### 4.1 use incar.job
        if opt_incar == 'j':
            if job and os.path.isfile(incar_job):
                incar = incar_job
            else:
                print(f"There is no {incar_job} in wdir")
        ### 4.2 use dir/INCAR
        elif os.path.isdir(opt_incar):
            incar = opt_incar + '/INCAR'
        ### 4.3 use directed INCAR file
        elif os.path.isfile(opt_incar):
            incar = opt_incar
    ### 4.4 iff INCAR is prepared in wdir
    elif job:
        if os.path.isfile(incar_job):
            incar = incar_job
    elif allprepared:
        print(f"INCAR in cwd will be used")
        incar = 'INCAR'
    ### INCAR in repo will not be used
    #elif os.path.isfile(incar_repo):
    #    incar = incar_repo 
    
    if "incar" not in locals():
        if os.path.isfile("INCAR"):
            incar = 'INCAR'
        else:
            print("Error:: cannot find INCAR")
            sys.exit()
    print(f'{incar} will be copied to INCAR')
    files2copy.append(f"{incar}")
    ### 5. make work_dir
    q = 'will you make dir? '
    if Lmkdir or yes_or_no(q):
        print(f"dirname {dirname}")
        if not "dirname" in locals():
            q = 'input dirname: '
            dirname = get_answers(q)
        if not os.path.isdir(dirname):
            com1 = "mkdir " + dirname
            print(com1)
            os.system(com1)
        else:
            q = f"{dirname} exists: want to overwrite?"
            if yes_or_no(q):
                pass
            else:
                print("Stop proceeding")
                sys.exit(1)
        
    for f in files2copy:
        if "INCAR" in f:
            com = f"cp {f} {dirname}/INCAR"
        else:
            com = f"cp {f} {dirname}"
        os.system(com)
    ### 6. check dir
    os.chdir(dirname)
    if not os.path.isfile('POTCAR'):
        if hostname == 'kisti':
            s = f"python {home}/sandboxg/pyvasp/genpotcar.py -pp pbe"
        else:
            #s = home + "/sandboxg/pyvasp/genpotcar.py -pp pbe"
            s = "genpotcar.py -pp pbe"
        if hpp_list:
            s += f" -hpp {' '.join(hpp_list)}"
        print(f"{s} in {dirname}")
        os.system(s)
    os.chdir(pwd)
    ##################################################
    ### run ? : first determin qx then qN
    if get_hostname()=='pt' and (not qx or not qN):
        qx, qN = get_queue_pt(qx=qx)
    s = qsub_command(dirname,X=qx,nnode=qN,np=qn, issue=issue)
    ### qsub runs inside module vas_qsub.qsub_command
    #print(f"{s}")
    #if Lrun or yes_or_no("Will you run ?"):
    #    os.system(s)

def main():
    parser = argparse.ArgumentParser(description='prepare vasp input files: -s for POSCAR -p POTCAR -k KPOINTS and -i INCAR')
    #parser.add_argument('-j', '--job', default="opt", choices=["pchg","chg","md","ini","zpe","mol","wav","opt","copt","sp", 'noD'], help='inquire for each file')
    parser.add_argument('-j', '--job', choices=["pchg","chg","md","ini","zpe","mol","wav","opt","copt","sp", 'noD'], help='inquire for each file')
    parser.add_argument('-s', '--poscar', help='poscar is required')
    parser.add_argument('-p', '--potcar', choices=['new','potpaw-pbe-new','old','potpaw-pbe-old','potpaw-gga'], help='pseudo potential directory: ')
    parser.add_argument('-hpp', '--pseudoH', nargs='*', help='include pseudo H list ')
    parser.add_argument('-k', '--kpoints', nargs='+', help='input number of k-points in kx, ky, kz, or g for gamma')
    ### toggle default: unset in the bare dir, set to j when INCAR.job exists
    parser.add_argument('-i', '--incar', help='j: use INCAR.job, dirname: use d/INCAR')
    parser.add_argument('-a', '--atoms', nargs='+', help='list of atoms')
    parser.add_argument('-f', '--iofile', default='incar.key', help='only read file is possible')
    parser.add_argument('-d', '--dname', help='get directory name')
    parser.add_argument('-al', '--all', action='store_true', help="prepared in job dir if not -s, -p, -k, -i")
    parser.add_argument('-r', '--run', action='store_true', help="submit job")
    parser.add_argument('-rd', '--mkdir', action='store_true', help="submit job")
    parser.add_argument('-err', '--error', choices=['opt','mem'], help="vasp error: converge, memory issue")
    g_queue = parser.add_argument_group(title='QUEUE')
    g_queue.add_argument('-x', '--xpartition', type=int, help="partition in platinum")
    g_queue.add_argument('-N', '--nnode', type=int, help="number of nodes, can be used to calculate total nproc")
    g_queue.add_argument('-n', '-np', '--nproc', help="number of nproc, total for pt, per node for kisti ")
    args = parser.parse_args()

    make_vasp_dir(args.job, args.poscar, args.potcar, args.pseudoH, args.kpoints, args.incar, args.all, args.dname, args.iofile, args.atoms, args.error, args.run, args.mkdir, args.xpartition, args.nnode, args.nproc)
    return 0

if __name__ == '__main__':
    main()
