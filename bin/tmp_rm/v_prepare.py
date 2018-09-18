#!/home/jackjack5/epd/bin/python

import argparse
import os
import glob
import re
import sys

def get_fkpoints(kname, sampling, kdir, n_kpoints):
    if kname:
        k_file  = kdir + '/kp.' + kname
    else:
        k_file  = kdir + '/kp.gamma'
    if not os.path.isfile(k_file):
        print ("there is not %s file" % (k_file))
        exit(31)
    cmd    = 'cp ' + k_file + ' ./KPOINTS'
    print ("KPOINTS was copied from %s to here" % k_file)
    return cmd

def get_fpp(pp, pp_dir, atoms):
    if atoms:
        cmd = 'cat '
        for atom in atoms:
            fname = pp_dir + '/' + atom +'.pot'
            if not os.path.isfile(fname):
                print("there is not %s.pot in %s" % (atom, potdir))
                exit(32)
            else:
                cmd += pp_dir + '/' + atom + '.pot '
        cmd += '  > ' + pp_dir + '/POTCAR'
        os.system(cmd)
        print ("POTCAR was generated from %s" % pp_dir)

    pot_file =  pp_dir + '/POTCAR'
    cmd1 = 'cp ' + pot_file + ' ./POTCAR'
    print ("POTCAR was copied from %s to here" % pot_file)
    return cmd1

def get_fincar(incar_name, in_path, incar_list):
    if not incar_name:
        incar_name='nonmag'
    in_file = in_path +'/incar.'+ incar_name
    if not os.path.isfile(in_file):
        print("there is not %s incar file" % (in_file))
        exit(33)
    cmd    = 'cp ' + in_file + ' ./INCAR'
    print ("INCAR was copied from %s to here" % in_file)
    return cmd

def get_jobs():
    ### obtain job_dir
    jobs = []
    jobs.extend(glob.glob('*.pos'))
    print jobs
    return jobs[0]

def modify_kpoints():
    pass
    return 0

def get_potcar(pp):
    pass
    return 0


def main():

    parser = argparse.ArgumentParser(description='prepare vasp input files: \n if -p: copy from VaspINI \n if not: use present directory')

    # job-kind: new or continous
    parser.add_argument('job', default='new', choices=['new','fcopy','cont','cont4','band','dos','pchg','soc'], help="kind of job: new or continuous \
    new: new job\
    fcopy: copy file without making a directory\
    others: copy from old job directory \
    ")

    # input job_dir
    parser.add_argument('job_dir', help='job directory name')

    # if new job
    parser.add_argument('--path', default='/qcfs/joonho/VaspINI', help='path to vasp initial input files, default=/qcfs/joonho/VaspINI/')
    
    # if cont job
    parser.add_argument('-o', '--odir', help='old dir for copying')
    parser.add_argument('--posd', help='POSCAR directory')

    #parser.add_argument('-f', '--sfiles', nargs='*', choices=['I','P','K','A'], help='define files to be copied from VaspINI dir')

    # what file to be copy?
    # INCAR
    groupi = parser.add_argument_group('incar')
    groupii=parser.add_mutually_exclusive_group()
    groupii.add_argument('-i', '--iname', default='nonmag',choices=['nonmag','vdw','vdw.so'], help='copy incar from default directory and by name')
    groupii.add_argument('--id', help='incar directory: default is Inc')
    groupi.add_argument('--functional', choices=['gga','ggapaw','pe','pepaw'], help='dft functional')
    # KPOINTS
    groupk = parser.add_argument_group('KPOINTS')
    groupkk= parser.add_mutually_exclusive_group()
    groupkk.add_argument('-k', '--kname', help='load KPOINTS using name, if not defined, gamma') 
    groupkk.add_argument('--kd', help='k-points directory, default is path/Kpoints')
    groupk.add_argument('--ksampling', default='mh', choices=['gamma', 'mh'], help='k-sampling method')
    groupk.add_argument('--knpoints', help='number of k-points')
    # POTCAR 
    groupp = parser.add_argument_group('pseudo potential')
    grouppp= parser.add_mutually_exclusive_group()
    grouppp.add_argument('-p', '--pname', default='potpaw_GGA', choices=['Pot-new', 'pot-old', 'potpaw_GGA'], help='pseudo potential directory')
    #grouppp.add_argument('--pname', choices=['gga', 'paw_pe', 'pe'], help='kind of pseudo-potential')
    grouppp.add_argument('-pd', help='pseudo potential file')
    groupp.add_argument('-a', '--atoms', nargs='+', help='list of atoms')
    args = parser.parse_args()
    
    in_default  = args.path + '/Inc'
    k_default   = args.path + '/Kpoints'
    p_default   = args.path 

    job         = args.job
    job_dir     = args.job_dir
    print args
    if not job == 'new' and not args.odir:
        print("Error & exit: input old directory using -o")
        print args.odir
        sys.exit(1)

    ### pwd::mk job directory
    if not job == 'fcopy':
        cmd = 'mkdir ' + job_dir
        if os.path.isdir(job_dir):
            print("%s directory is already there" % (job_dir))
            exit(10)       # should exit when being used, now testing
        else:        
            os.system(cmd)

    odir_files=['POSCAR', 'KPOINTS', 'POTCAR', 'INCAR', 'WAVECAR', 'CHGCAR']

    incar_list=['start','mag','prec','parallel', 'gga', 'scf']
    # start vs cont
    # mag vs nonmag
    # scf vs post-scf
    # post-scf: dos, band, pchg
    # soc ?

    cmds=[]

    ### if new job
    if job == "new":
        ### 1 pwd::copy POSCAR        
        st_poscar = job_dir + '.pos'        
        if os.path.isfile(st_poscar):
            cmd = 'cp ' + st_poscar + ' ' + job_dir + '/POSCAR'
            os.system(cmd)
            print ("%s was copied to %s" % (st_poscar, job_dir))
        else:
            print ("Warning::POSCAR cannot be found by %s" % (st_poscar))

        ### path check for initial upload of files            
        if not os.path.isdir(args.path):
            print ("there is no path to %s" % (args.path))
            exit(30)
        ### path2vasp_ini:: 2 copy KPOINTS
        if args.kd:
            kpath = args.kd
        else:
            kpath = args.path + '/Kpoints'
        cmd1 = get_fkpoints(args.kname, args.ksampling, kpath, args.knpoints)
        cmds.append(cmd1)
        cmd1    = 'cp KPOINTS ' + job_dir
        cmds.append(cmd1)
        ### path2vasp_ini:: 3 copy pseudo-potential
        pp_path = args.path + '/' + args.pd
        if not os.path.isdir(pp_path):
            print("there is not dir for POTCAR in %s" % (pp_path))
            exit(320)
        cmd2 = get_fpp(args.pname, pp_path, args.atoms)
        cmds.append(cmd2)
        cmd2    = 'cp POTCAR ' + job_dir
        cmds.append(cmd2)
        ### path2vasp_ini:: 4 copy INCAR
        in_path = args.path +'/'+ args.idir
        cmd3 = get_fincar(args.incar, in_path, incar_list)
        cmds.append(cmd3)
        cmd3    = 'cp INCAR ' + job_dir
        cmds.append(cmd3)
    elif job == 'fcopy':
        #for file in args.sfiles:
        if args.kname:
            kpath = args.path + '/Kpoints'
            cmd1 = get_fkpoints(args.kname, args.ksampling, kpath, args.knpoints)
            cmds.append(cmd1)
            cmd1    = 'cp KPOINTS ' + job_dir
            cmds.append(cmd1)
        if args.pname:
            cmd2 = get_fpp(args.pname, pp_path, args.atoms)
            cmds.append(cmd2)
            cmd2    = 'cp POTCAR ' + job_dir
            cmds.append(cmd2)
        if args.incar:
            cmd3 = get_fincar(args.incar, in_path, incar_list)
            cmds.append(cmd3)
            cmd3    = 'cp INCAR ' + job_dir
            cmds.append(cmd3)

    ### if continue, copy from old-directory
    else:
        if not args.odir:
            print ('job== %s requires odir' % (args.odir))
            exit(2)
        print ("copy from %s to  %s" % (args.odir, job_dir))

        if job == 'cont4':
            files=odir_files[:4]
        else:
            files=odir_files[:]
        for cpfile in files:
            if cpfile == 'POSCAR':
                filename = args.odir + '/CONTCAR'
                if os.path.isfile(filename):
                    old_file = 'CONTCAR'
                else:
                    old_file = 'POSCAR'
            elif cpfile == 'KPOINTS':
                filename = args.odir + '/IBZKPT'
                if os.path.isfile(filename):
                    old_file = 'IBZKPT'
                else:
                    old_file = 'KPOINTS'
            else:
                old_file = cpfile
            filename = args.odir + '/' + old_file
            cmd0    = 'cp ' + filename + '  ' + job_dir + '/' + cpfile
            os.system(cmd0)
            print ("%s is copied from %s to %s" % (old_file, args.odir, job_dir))

        if args.kname:
            kpath = args.path + '/Kpoints'
            cmd = get_fkpoints(args.kname, args.ksampling, kpath, args.knpoints)
            os.system(cmd)
            cmd1    = 'cp KPOINTS ' + job_dir
            cmds.append(cmd1)

    ### if id, kd, pd is defined, load from the input directory and overwrite
    if args.id:
        cmdn = 'cp '+args.id+ '/INCAR '+job_dir 
        cmds.append(cmdn)
        print ("INCAR is copied from %s to %s" % (args.id, job_dir))
    if args.kd:
        filename = args.kd + '/IBZKPT'
        if os.path.isfile(filename):
            f = "IBZKPT"
        else:
            f = "KPOINTS"
        cmdn = 'cp '+ args.kd + '/' + f
        cmds.append(cmdn)
        print ("%s is copied from %s to %s" % (f, args.kd, job_dir))
    if args.pd:
        cmdn = 'cp '+args.pd+ '/POTCAR ' + job_dir
        cmds.append(cmdn)
        print ("POTCAR is copied from %s to %s" % (args.pd, job_dir))
    if args.posd:
        cmdn = 'cp '+args.posd+'/POSCAR ' + job_dir
        cmds.append(cmdn)
        print ("POSCAR is copied from %s to %s" % (args.posd, job_dir))

    for command in cmds:            
        os.system(command)
            
    return 0


if __name__ == '__main__':
    main()
