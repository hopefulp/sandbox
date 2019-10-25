#!/home/jackjack5/epd/bin/python

import argparse
import os
import glob
import re

def get_jobs():
    ### obtain new_dir
    jobs = []
    jobs.extend(glob.glob('*.pos'))
    print(jobs)
    return jobs[0]

def modify_kpoints():
    pass
    return 0

def get_potcar(pp):
    pass
    return 0


def main():

    parser = argparse.ArgumentParser(description='prepare vasp input files: \n if -p: copy from VaspINI \n if not: use present directory')
    parser.add_argument('job', default='new', choices=['new','cont', 'cont4','band'], help='what kind of directory')
    parser.add_argument('new_dir', help='new directory name')

    groupj = parser.add_argument_group('jobs')
    groupj.add_argument('-f', '--files', nargs='*', choices=['I','P','K','A'], default='N', help='define files to be copied from VaspINI dir')

    groupo = parser.add_mutually_exclusive_group()
    groupo.add_argument('-p', '--path', action='store_const', const='/qcfs/joonho/VaspINI/', help='path to vasp initial input files, default=/qcfs/joonho/VaspINI/; if not evoked, files in pwd will be copied to job directory')
    groupo.add_argument('-o', '--odir', help='dir for copy initial files')
    
    groupi = parser.add_argument_group('incar')
    groupi.add_argument('-i', '--incar', default='nonmag', choices=['nonmag','vdw','vdw.so'], help='copy incar')
    groupi.add_argument('--idir', default='Inc', choices=['Inc', 'INC'], help='incar directory')
    groupi.add_argument('--functional', choices=['gga','ggapaw','pe','pepaw'], help='dft functional')


    groupk = parser.add_argument_group('k-points')
    groupk.add_argument('-k', '--ksample', default='mh', choices=['gamma', 'mh'], help='k-sampling method')
    groupk.add_argument('--kdir', default='Kpoints', help='k-points directory')
    groupk.add_argument('--kpoints', help='number of k-points')
    

    groupp = parser.add_argument_group('pseudo potential')
    groupp.add_argument( '--pp', choices=['gga', 'paw_pe', 'pe'], help='kind of pseudo-potential')
    groupp.add_argument('--pdir', default='potpaw_GGA', choices=['Pot-new', 'pot-old', 'potpaw_GGA'], help='pseudo potential directory')
    groupp.add_argument('-a', '--atoms', nargs='+', help='list of atoms')
    args = parser.parse_args()

    cp_files=['CONTCAR', 'IBZKPT', 'POTCAR', 'INCAR', 'WAVECAR', 'CHGCAR']
    if re.search('cont', args.job):
        if not args.odir:
            print(('job== %s requires odir' % (args.odir)))
            exit(2)
        print(("copy from %s to  %s" % (args.odir, args.new_dir)))
        cmd0 = 'mkdir ' + args.new_dir
        os.system(cmd0)
        if args.job == 'cont':
            files=flist_all[:]
        elif args.job == 'cont4':
            files=flist_all[:4]
        for cfile in files:
            if cfile == 'CONTCAR':
                obj =  'POSCAR'
            else:
                obj = ''
            filename = args.odir + '/' + cfile
            if not os.path.isfile(filename):
                print(('there is not %s' % (filename)))
                exit(3)
            cmd0    = 'cp ' + filename + '  ' + args.new_dir + '/' + obj
            os.system(cmd0)
        return 0                   

    if args.new_dir:
        new_dir = args.new_dir
    else:        
        # only one job is assigned
        poscar  = get_jobs()
        job     = []
        job     = poscar.split('.')
        new_dir = job[0]
        #print new_dir
            
    ### pwd::mk job directory
    cmd = 'mkdir ' + new_dir
    if os.path.isdir(new_dir):
        print(("%s directory is already there" % (new_dir)))
        exit(10)       # should exit when being used, now testing
    else:        
        os.system(cmd)
    ### 1 pwd::copy POSCAR        
    st_poscar = new_dir + '.pos'        
    if os.path.isfile(st_poscar):
        cmd = 'cp ' + st_poscar + ' ' + new_dir + '/POSCAR'
        os.system(cmd)
        print(("%s was copied to %s" % (st_poscar, new_dir)))
    #else:
    #    cmd = 'cp POSCAR ' + new_dir
    #    os.system(cmd)
    #    print("POSCAR was copied to %s" % (new_dir))

    ### COPY files
    ### if -p is evoken, use VASP-INI
    #if args.files == 'N':
    #    return 0
    cmds=[]
    if args.path:

        if not os.path.isdir(args.path):
            print(("there is no path to %s" % (args.path)))
            exit(30)
        ### path2vasp_ini:: 2 copy KPOINTS
        if 'K' in args.files or 'A' in args.files:

            k_file  = args.path + '/kp.' + args.ksample
            if not os.path.isfile(k_file):
                print(("there is not %s file" % (k_file)))
                exit(31)
            cmd1    = 'cp ' + k_file + ' ./KPOINTS'
            cmds.append(cmd1)
            print("KPOINTS was copied from VaspINI to here")
            cmd1    = 'cp KPOINTS ' + new_dir
            cmds.append(cmd1)
        ### path2vasp_ini:: 3 copy pseudo-potential
        if 'P' in args.files or 'A' in args.files:
            pot_dir = args.path +'/'+ args.pdir
            if not os.path.isdir(pot_dir):
                print(("there is not dir for POTCAR in %s" % (args.pdir)))
                exit(320)
            # make POTCAR                
            if args.atoms:
                cmd = 'cat '
                for at in args.atoms:
                    fname = pot_dir + '/' + at +'.pot'
                    if not os.path.isfile(fname):
                        print(("there is not %s.pot in %s" % (at, potdir)))
                        exit(32)
                    else:
                        cmd += pot_dir + '/' + at + '.pot '
                cmd += '  > ' + pot_dir + '/POTCAR'
                os.system(cmd)
            pot_file =  pot_dir + '/POTCAR'
            #if not os.path.isfile(pot_file):
            #    print("there is not POTCAR in %s" % (args.pdir))
            #    exit(32)
            cmd2    = 'cp ' + pot_file + ' ./POTCAR'
            cmds.append(cmd2)
            print("POTCAR was copied from VaspINI to here")
            cmd2    = 'cp POTCAR ' + new_dir
            cmds.append(cmd2)
        if 'I' in args.files or 'A' in args.files:            
            in_file = args.path +'/'+ args.idir +'/incar.'+ args.incar
            if not os.path.isfile(in_file):
                print(("there is not %s incar file" % (in_file)))
                exit(33)
            cmd3    = 'cp ' + in_file + ' ./INCAR'
            cmds.append(cmd3)
            print("INCAR was copied from VaspINI to here")
            cmd3    = 'cp INCAR ' + new_dir
            cmds.append(cmd3)
    ### if not -p: copy files to working-dir 
    elif args.odir:
        print(("copy from %s" % (args.odir)))
        cmd1    = 'cp ' + args.odir + '/KPOINTS ' + new_dir
        cmd2    = 'cp ' + args.odir + '/POTCAR ' + new_dir
        cmd3    = 'cp ' + args.odir + '/INCAR ' + new_dir
        cmds.append(cmd1)
        cmds.append(cmd2)
        cmds.append(cmd3)
        file    = new_dir + '/POSCAR'
        if not os.path.isfile(file):
            cmd4    = 'cp ' + args.odir + '/POSCAR ' + new_dir
            cmds.append(cmd4)
    else:            
        ### 2 copy KPOINTS
        if os.path.isfile('KPOINTS'):
            cmd1    = 'cp ' + 'KPOINTS ' + new_dir
            cmds.append(cmd1)
        else:
            print ("there is not KPOINTS")
            exit(41)

        ### 3 copy pseudo-potential
        if os.path.isfile('POTCAR'):
            cmd2    = 'cp ' + 'POTCAR ' + new_dir
            cmds.append(cmd2)
        else:
            print ("there is not POTCAR")
            exit(42)
        ### path2vasp_ini:: 4 copy INCAR
        if os.path.isfile('INCAR'):
            cmd3    = 'cp ' + 'INCAR ' + new_dir
            cmds.append(cmd3)
        else:
            print ("there is not INCAR")
            exit(43)
    for command in cmds:            
        os.system(command)
            
    ### modify k-points
    if args.kpoints:
        modify_kpoints()

    return 0


if __name__ == '__main__':
    main()
