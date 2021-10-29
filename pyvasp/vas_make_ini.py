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
from  myvasp import *
from common import *
from vas_qsub import qsub_command
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

def make_vasp_dir(job, poscar, apotcar, kpoints, incar, allprepared, Lquestion, kpsub, dirname, iofile, atoms, Lrun, qopt):
    global ini_dvasp, pwd
    ### 0. obtain default vasp repository
    ini_dvasp = get_vasp_repository()
    pwd = os.getcwd()
    files2copy=[]
    ### Now this is running
    if Lquestion:
        ### 1. get POSCAR: make dirname using poscar
        if not poscar:
            q = 'will you make POSCAR? '
            if yes_or_no(q):
                q = 'input file: '
                poscar = get_answers(q)
        else:
            ### cp input poscar to 'POSCAR'
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
        if allprepared:
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
        files2copy.append('KPOINTS')
        ### 4. get INCAR :: use make_incar.py
        q = 'will you make INCAR? '
        incarjob = 'INCAR.' + job
        q2 = f'will you use {incarjob}?'
        if allprepared:
            print(f"INCAR in cwd will be copied to {dirname}")
            incar = 'INCAR'
        elif os.path.isfile(incarjob) and yes_or_no(q2):
            incar = incarjob
            os.system(f'cp {incar} INCAR')
        elif yes_or_no(q):
            q = 'input incar-key file or use make_incar.py: '
            keyfile = get_answers(q)
            if not keyfile:
                keyfile = iofile
            get_incar(keyfile)
        files2copy.append('INCAR')
        ### 5. make work_dir
        q = 'will you make dir? '
        if yes_or_no(q):
            print(f"dirname {dirname}")
            if not "dirname" in locals():
                q = 'input dirname: '
                dirname = get_answers(q)
            if not os.path.isdir(dirname):
                com1 = "mkdir " + dirname
                print(com1)
                os.system(com1)
            
        for f in files2copy:
            com2 = f"cp {f} " + dirname
            os.system(com2)
        ### 6. check dir
        os.chdir(dirname)
        if not os.path.isfile('POTCAR'):
            if hostname == 'kisti':
                s = f"python {home}/sandboxg/pyvasp/genpotcar.py -pp pbe"
            else:
                #s = home + "/sandboxg/pyvasp/genpotcar.py -pp pbe"
                s = "genpotcar.py -pp pbe"
            os.system(s)
            print(f"in {dirname}")
        os.chdir(pwd)
    ### without -q            
    else:
        ### 1. get POSCAR
        dirname = 'tmp_vasp'
        if aposcar:
            get_poscar(aposcar)
        else:
            get_poscar("")

        ### 2. get POTCAR
        if apotcar:
            if not atoms:
                q = 'input atoms in the order of poscar: '
                atoms = get_answers(q).split()
                get_potcar(apotcar, atoms)
            else:            
                get_potcar(apotcar, atoms)
        else:
            if not atoms:
                print("Use -a atom list as minimum requirement")
            else:
                get_potcar('new', atoms) 
            
        ### 3. get KPOINTS
        if kpoints:
            make_kpoints(kpoints, kpsub)
        else:
            print("KPOINTS will be made from gamma")
            make_kpoints("", 'gamma')
        ### 4. get INCAR :: use make_incar.py        
        if iofile:
            print("Used 'incar.key'")
            get_incar("incar.key")
    ##################################################
    ### run ?
    if Lrun:
        s = qsub_command(dirname, qopt=qopt)
        print(f"{s}")
        os.system(s)
        


def main():
    parser = argparse.ArgumentParser(description='prepare vasp input files: -s for POSCAR -p POTCAR -k KPOINTS and -i INCAR')
    parser.add_argument('-j', '--job', default="opt", choices=["dos","band","pchg","chg","md","ini","zpe","mol","wav","opt","sp"], help='inquire for each file')
    parser.add_argument('-q', '--question', action='store_false', help='inquire for each file')
    parser.add_argument('-s', '--poscar', help='poscar is required')
    parser.add_argument('-p', '--potcar', choices=['new','potpaw-pbe-new','old','potpaw-pbe-old','potpaw-gga'], help='pseudo potential directory: ')
    parser.add_argument('-a', '--atoms', nargs='+', help='list of atoms')
    parser.add_argument('-k', '--kpoints', nargs='+', help='input number of k-points in kx, ky, kz')
    parser.add_argument('-ks', '--kpsub', default='monk', choices=['monk','gamma','dos','band'], help='diverse k-point sampling')

    parser.add_argument('-i', '--incar', action='store_true',  help='first run make_incar.py then use incar.key')
    parser.add_argument('-f', '--iofile', default='incar.key', help='only read file is possible')
    parser.add_argument('-d', '--directory', help='mkdir and cp')
    parser.add_argument('-al', '--all', action='store_true', help="skip if KPOINTS, POTCAR, INCAR were ready")
    parser.add_argument('-o', '--qopt', action='store_true', help="use different qsub file in kisti")
    parser.add_argument('-r', '--run', action='store_true', help="submit job")
    args = parser.parse_args()

    make_vasp_dir(args.job, args.poscar, args.potcar, args.kpoints, args.incar, args.all, args.question, args.kpsub, args.directory, args.iofile, args.atoms, args.run, args.qopt)
    return 0

if __name__ == '__main__':
    main()
