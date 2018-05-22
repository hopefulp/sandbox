#!/usr/bin/python
### versin 1.1 by J. Park
### 2018.4.2 makes input files by option -s(POSCAR) -p(POTCAR) -k(KPOINTS) -i(INCAR)
### incar is not ready

import argparse
import os
import shutil
import re
from  myvasp import *
from common import *

pseudo_pot={'new':'Pot-new', 'potpaw-pbe-new':'Pot-new', 'old':'pot-old', 'potpaw-pbe-old':'pot-old'}
#global pwd, ini_dvasp

def get_poscar(poscar):
    # confirm file location
    if not os.access('%s' % poscar, os.F_OK):
        print 'poscar is not detectable'
        exit(2)
    else:
        comm = 'cp %s POSCAR' % poscar
        os.system(comm)
        print 'POSCAR is made'
    # confirm POSCAR is made        
    if not os.access('POSCAR', os.F_OK):
        print 'POSCAR is not here'
        exit(21)
            
    return 0    
def get_answers(question):
    reply = str(raw_input(question)).strip()
    return reply

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
            print "atom %s is not detected" % atom
            exit(32)
    comm = 'cat '            
    for potfile in potfiles:
        fname = potdir + '/' + potfile
        comm += fname + ' '
    comm += ' > POTCAR'
    os.system(comm)
    print 'POTCAR is combined from atoms ', atoms

    return 0    

def get_incar(ifile):
    dic = {}
    rw = 'r'
    make_incar(dic, 'r', ifile)
    print 'INCAR is made from %s' % ifile

    return 0

def main():
    global ini_dvasp, pwd
    parser = argparse.ArgumentParser(description='prepare vasp input files: -s for POSCAR -p POTCAR -k KPOINTS and -i INCAR')
    parser.add_argument('-q', '--question', action='store_true', help='inquire for each file')
    parser.add_argument('-s', '--poscar', help='poscar is required')
    parser.add_argument('-p', '--potcar', choices=['new','potpaw-pbe-new','old','potpaw-pbe-old','potpaw-gga'], help='pseudo potential directory: ')
    parser.add_argument('-a', '--atoms', nargs='+', help='list of atoms')
    parser.add_argument('-k', '--kpoints', nargs='+', help='input number of k-points in kx, ky, kz')
    parser.add_argument('-l', '--ksub', default='monk', choices=['monk','gamma','dos','band'], help='diverse k-point sampling')

    parser.add_argument('-i', '--incar', action='store_true',  help='first run make_incar.py then use incar.key')
    parser.add_argument('-f', '--iofile', default='incar.key', help='only read file is possible')
    parser.add_argument('-d', '--directory', help='mkdir and cp')
    args = parser.parse_args()

    ### 1. obtain default vasp repository
    ini_dvasp = get_vasp_repository()
    pwd = os.getcwd()
    if args.question:
        ### 2. get POSCAR
        q = 'will you make POSCAR? '
        if yes_or_no(q):
            q = 'input .pos file: '
            poscar = get_answers(q)
            get_poscar(poscar)
        ### 3. get POTCAR
        q = 'will you make POTCAR? '
        if yes_or_no(q):
            q = 'input pseudo-potential type (new, old, gga, etc): '
            pot = get_answers(q)
            q = 'input atoms in the order of poscar: '
            atoms = get_answers(q).split()
            get_potcar(pot, atoms)
        ### 4. get KPOINTS
        q = 'will you make KPOINTS?'
        if yes_or_no(q):
            q = 'input nummber of kpoints: '
            kp = get_answers(q).split()
            if len(kp) == 4:
                method = kp.pop(3)
            elif len(kp) == 3:
                method = 'monk'
                print 'default is MH'
            print kp, method
            make_kpoints(kp, method)
        ### 5. get INCAR :: use make_incar.py
        q = 'will you make INCAR? '
        if yes_or_no(q):
            q = 'give input incar-key file or use make_incar.py: '
            keyfile = get_answers(q)
            get_incar(keyfile)
            
    else:
        ### 2. get POSCAR
        if args.poscar:
            get_poscar(args.poscar)
        ### 3. get POTCAR
        if args.potcar:
            if not args.atoms:
                q = 'input atoms in the order of poscar: '
                atoms = get_answers(q).split()
                get_potcar(args.potcar, atoms)
            else:            
                get_potcar(args.potcar, args.atoms)
        ### 4. get KPOINTS
        if args.kpoints:
            make_kpoints(args.kpoints, args.ksub)
        ### 5. get INCAR :: use make_incar.py        
        if args.incar:
            get_incar(args.iofile)
                
    ### 6. mkdir and cp POSCAR POTCAR KPOINTS INCAR
    if args.directory:
        os.mkdir(args.directory)
        print "directory ./%s was made" % args.directory
        shutil.copy('POSCAR', args.directory)
        shutil.copy('POTCAR', args.directory)
        shutil.copy('KPOINTS', args.directory)
        shutil.copy('INCAR', args.directory)
        print 'POSCAR POTCAR KPOINTS INCAR were copied'
        

if __name__ == '__main__':
    main()
