#!/usr/bin/python

import argparse
import os
import subprocess


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
    
    print "vasp repository is ", ini_dvasp, ' in system ', hostname
    if not os.access(ini_dvasp, os.F_OK):
        print "Error:: the directory cannot be found\n stop"
        exit(1)
    return 0

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

def get_potcar(pot,atoms):
    potfiles=[]
    if pot == 'None':
        if not os.access('./POTCAR', os.F_OK):
            print 'POTCAR is not here'
            exit(3)
        return 0
    if not atoms:
        print 'natoms are required'
        exit(31)
    potdir = ini_dvasp + '/' + pot
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
            


def main():
    parser = argparse.ArgumentParser(description='prepare vasp input files: try poscar and atom list')
    parser.add_argument('poscar', help='poscar is required')
    parser.add_argument('atoms', nargs='+', help='list of atoms')
    parser.add_argument('-p', '--potcar', default='potpaw-pbe-new', choices=['potpaw-pbe-new', 'potpaw-pbe-old', 'potpaw-gga', 'None'], help='pseudo potential directory: ')
    args = parser.parse_args()

    ### 1. obtain default vasp repository
    get_vasp_repository()

    ### 2. get POSCAR
    get_poscar(args.poscar)

    ### 3. get POTCAR
    get_potcar(args.potcar, args.atoms)

    ### 4. get KPOINTS
    get_kpoints()    

if __name__ == '__main__':
    main()
