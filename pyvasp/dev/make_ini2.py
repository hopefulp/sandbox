#!/usr/bin/python
### versin 1.0 by J. Park

import argparse
import os
import subprocess

pseudo_pot={'new':'Pot-new', 'potpaw-pbe-new':'Pot-new', 'old':'pot-old', 'potpaw-pbe-old':'pot-old'}


def get_vasp_repository():

    global ini_dvasp
    hostname = os.popen('hostname').read().rstrip()
    #print 'hostname is', hostname

    ### subprocess is not working
    #hostname = subprocess.popen('hostname', stdout=subprocess.PIPE, shell=True)
    #proc = subprocess.Popen(["cat", "/etc/services"], stdout=subprocess.PIPE, shell=True)
    #(out, err) = proc.communicate()
    #print "program output:", out
    if hostname == 'chi':
        ini_dvasp = '/Data/Bkaist/VaspINI'
    else:
        ini_dvasp = '/qcfs/joonho/VaspINI'
    
    print("vasp repository is ", ini_dvasp, ' in system ', hostname)
    if not os.access(ini_dvasp, os.F_OK):
        print("Error:: the directory cannot be found\n stop")
        exit(1)
    return 0

def yes_or_no(question):
    reply = str(input(question+' (y/n): ')).lower().strip()
    if reply[0] == 'y':
        return True
    else:
        print('skipped')
        return False

def get_poscar(poscar):
    # confirm file location
    if not os.access('%s' % poscar, os.F_OK):
        print('poscar is not detectable')
        exit(2)
    else:
        comm = 'cp %s POSCAR' % poscar
        os.system(comm)
        print('POSCAR is made')
    # confirm POSCAR is made        
    if not os.access('POSCAR', os.F_OK):
        print('POSCAR is not here')
        exit(21)
            
    return 0    
def get_atoms(question):
    reply = str(input(question)).split()
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

def get_kpoints(kp,method):
    kpdir = ini_dvasp + '/' + 'Kpoints/'
    fname = 'kp'+ kp +'.'+method
    fullname = kpdir+fname
    if os.access(fullname, os.F_OK):
        comm = 'cp '+ fullname + ' ./KPOINTS'
        os.system(comm)
        print('KPONTS was copied')
    else:
        print(fullname, ' is not detectable')
        exit(41)
    return 0
     
def get_incar(scf, nonscf):
    incars=[]
    idir = ini_dvasp + '/' + 'INC/'
    
    return 0


def main():
    parser = argparse.ArgumentParser(description='prepare vasp input files: -s for POSCAR -p POTCAR -k KPOINTS and -i INCAR')
    parser.add_argument('-s', '--poscar', help='poscar is required')
    parser.add_argument('-p', '--potcar', choices=['new','potpaw-pbe-new','old','potpaw-pbe-old','potpaw-gga'], help='pseudo potential directory: ')
    parser.add_argument('-a', '--atoms', nargs='+', help='list of atoms')
    parser.add_argument('-k', '--kpoints', help='input number of k-points in kx, ky, kz')
    parser.add_argument('-l', '--ksub', default='monk', choices=['monk','gamma','dos','band'], help='diverse k-point sampling')

    parser.add_argument('-i', '--incar', default='+', choices=['ini','cont','mag','afm','pbe'], help='diverse incar contents')
    parser.add_argument('-j', '--nonscf', default='+', choices=['dos','band','pchg','soc'],help='incar for non-scf calculation')
    args = parser.parse_args()

    ### 1. obtain default vasp repository
    get_vasp_repository()
    ### 2. get POSCAR
    if args.poscar:
        get_poscar(args.poscar)
    ### 3. get POTCAR
    if args.potcar:
        if not args.atoms:
            atoms = get_atoms("input atoms in the order of poscar: ")
            get_potcar(args.potcar, atoms)
        else:            
            get_potcar(args.potcar, args.atoms)
    ### 4. get KPOINTS
    if args.kpoints:
        get_kpoints(args.kpoints, args.ksub)
    if args.incar:
        get_incar(args.incar, args.nonscf)

if __name__ == '__main__':
    main()
