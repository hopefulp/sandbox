#!/home/joonho/anaconda3/bin/python
### versin 1.1 by J. Park
### 2018.4.2 makes input files by option -s(POSCAR) -p(POTCAR) -k(KPOINTS) -i(INCAR)
### incar is not ready
### 2019.10.25 update
### 2021.05.07 update for EE
### 2023.04.11 poscar canbe put by prefix: if not make dir, poscar is skipped

import argparse
import os, sys
import shutil
import re
import string
from common     import *
from vas_qsub   import get_queue_pt, qsub_command
from mod_vas    import *
from libposcar import get_poscar, get_dnames4pos 
from libincar  import modify_incar_bykv, add_inckv_bysubjob

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
################# 1       2        3       4         5        6     7           8         9        10      11    12 13 14   15   16  
def make_vasp_dir(job, subjob, poscars, apotcar, jobadds, kpoints, incar, incopt, dirnames, option, allprepared, iofile, qx,qN,qn,vasp_exe,lkisti,Lrun):
    '''
    job
    subjob  in case of md, nve, NVT: heat, cool 
    poscars     list of poscar
    apotcar
    jobadds NEB final POSCAR
            ZPE atom lists
            ... pseudo potentials for H
    kpoints     list of 3 digits (str)
    incar
    dirnames
    option
    allprepared
    '''
    global ini_dvasp, pwd

    if job == 'neb':
        poscar_fin = jobadds[0]
        if len(jobadds) == 2:
            nimages = jobadds[1]
    elif job == 'pseudo':
        hpp_list = jobadds
    elif job == 'zpe':
        atoms = jobadds

    ### 0. obtain default vasp repository
    ini_dvasp = get_vasp_repository()
    pwd = os.getcwd()
    #if not os.path.isfile(poscars[0]):
    #    print(f"can't find POSCAR")
    #    sys.exit(1)
    print(f"{poscars} {dirnames}")
    for poscar, dirname in zip(poscars, dirnames):
        print(f"target directory {dirname}")
        ### 1.1 make work_dir
        #q = f'will you make dir? {dirname}'
        #if yes_or_no(q):
        ### kpoins are input one by one
        #if Lktest:
        #    dirname += 'k' + list2str(kpoints)
        #print(f"dirname {dirname} in function {whereami()}()")
        if not os.path.isdir(dirname):
            com1 = "mkdir " + dirname
            print(com1)
            os.system(com1)
        #### if dir exists
        else:
            q = f"{dirname} exists: want to overwrite?"
            #if Lrun or  yes_or_no(q):
            ### overwrite
            if Lrun and ('a' in Lrun or 'o' in Lrun):
                pass
            ### No modification & sumbit and go to next ele in for-loop
            elif Lrun and 's' in Lrun:
                s = qsub_command(dirname,X=qx,nnode=qN,np=qn, issue=issue, vasp_exe=vasp_exe, lkisti=lkisti, Lrun=Lrun)
                print(f"Job {dirname}  was submitted without modification")
                continue
            ### overwrite or not
            elif yes_or_no(q):
                pass
            ### stop for this directory
            else:
                print("skip overwriting")
                print("Go to next poscar")
                continue
        ### overwrite continues
        ### 1.2 Copy POSCAR
        try :
            os.path.isfile(poscar)
            com = f"cp {poscar} {dirname}/POSCAR"
            print(f"{com}")
            os.system(f'{com}')
        ### this except is not working: OSError
        except FileNotFoundError :
            print(f"cannot find {poscar}")
            sys.exit(11)
        #else:
        #    print("POSCAR error")
        #    sys.exit(12)
            
        
        ### 2. get KPOINTS
        ### use kp inputs, if KPOINTS ready, KPOINTS.job, make KPOINTS file
        q = 'will you make KPOINTS?'

        if job and os.path.isfile(f"KPOINTS.{job}"):
            com = f'cp KPOINTS.{job} {dirname}/KPOINTS'
            os.system(f'{com}')
            print(f"{com}")
        elif allprepared:
            com = f'cp KPOINTS {dirname}'
            os.system(f'{com}')
            print(f"KPOINTS in cwd was copied to {dirname}")
        elif kpoints:
            if len(kpoints) == 3 and kpoints[0] == '1' and kpoints[1] == '1' and kpoints[2] =='1':
                kpoints = []
                kpoints.append('g')
            if len(kpoints) == 3: 
                method = "MH"
            if len(kpoints) == 1 and kpoints[0] == 'g':
                kpoints = "1  1  1"
                method = "gamma"
            make_kpoints(kpoints, method)
            com = f'cp KPOINTS {dirname}'
            os.system(f'{com}')
            print(f"KPOINTS was made and copied to {dirname}")
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
            com = f'cp KPOINTS {dirname}'
            os.system(f'{com}')
            print(f"KPOINTS was made and copied to {dirname}")
        else:
            com = f'cp KPOINTS {dirname}'
            os.system(f'{com}')
            print("Use wdir KPOINTS")

        ### 3. get INCAR :: use make_incar.py
        incar_repo=f"{ini_dvasp}/INCAR.{job}"
        ### 3.1 Use inserted INCAR by -i option
        
        if incar and os.path.exists(incar):
            if os.path.isfile(incar):
                f_incar = incar
            elif os.path.isdir(incar):
                f_incar = incar + '/INCAR'
        ### 3.2 if job exists, use INCAR.job
        elif job and os.path.isfile(f"INCAR.{job}"):
            f_incar = f"INCAR.{job}"
        ### 3.3 if INCAR is prepared in wdir
        elif allprepared and os.path.isfile("INCAR"):
            print(f"INCAR in cwd will be used")
            f_incar = 'INCAR'
        else:
            print("Error:: cannot find INCAR")
            sys.exit()

        ### Modify INCAR
        ### subjob for MD: cool, heat, quench
        if job == 'md' and subjob:
            print(f"Add more INCAR KW for {subjob} if not specified in incopt {incopt}")
            if not incopt:
                incopt = {}
            else:
                incopt = list2dict(incopt)
            add_inckv_bysubjob(job, subjob, incopt)
            print(f"incopt after add_inckv_bysubjob {incopt}")
        ### though job != 'md', incopt would work
        if incopt:
            if type(incopt) != dict: 
                print(f"1st change of incopt {incopt}")
                incopt = list2dict(incopt)
            modify_incar_bykv(f_incar, incopt, outf='INCAR.new', mode='m')
            f_incar = 'INCAR.new'

        com = f'cp {f_incar} {dirname}/INCAR'
        os.system(f'{com}')
        print(f"{f_incar} was copied to {dirname}/INCAR")

        ###### 6. make POTCAR; job=neb: cp POSCAR_fin, extract NIMAGES
        ### {wdir} NEB: cp POSCAR_fin
        if job == 'neb':
            if 'poscar_fin' in locals():
                s = f'cp {poscar_fin} {dirname}/POSCAR_fin'
                os.system(s)
                print(f"{poscar_fin} is added for final POSCAR")
            else:
                print(f'final POSCAR is required by -ja POSCAR_fin')
                sys.exit(99)

        ### {jobdir} NEB: extract nimages from INCAR; Make POTCAR
        os.chdir(dirname)
        if job == 'neb':
            nimages = modify_incar_bykv('INCAR', ['IMAGES'], mode='e')[0]
            s = f'nebmake.pl POSCAR POSCAR_fin {nimages}'
            os.system(s)
        #if not os.path.isfile('POTCAR'):   # to make new POTCAR 
        if hostname == 'kisti':
            s = f"python {home}/sandboxg/pyvasp/genpotcar.py -pp pbe"
        else:
            #s = home + "/sandboxg/pyvasp/genpotcar.py -pp pbe"
            s = "genpotcar.py -pp pbe"
        if 'hpp_list' in locals():
            s += f" -hpp {' '.join(hpp_list)}"
        print(f"{s} in {dirname}")
        os.system(s)
        os.chdir(pwd)
        ##################################################
        ### first determine qx then qN for pt
        ### if Lrun == n, pass but None, ask in qsub_command
        if Lrun and 'n' in Lrun:
            pass
        else:
            if get_hostname()=='pt' and (not qx or not qN):
                qx, qN = get_queue_pt(qx=qx)
            s = qsub_command(dirname,X=qx,nnode=qN,np=qn, option=option, vasp_exe=vasp_exe, lkisti=lkisti, Lrun=Lrun)

    return 0            

def main():
    parser = argparse.ArgumentParser(description='prepare vasp input files: -s for POSCAR -p POTCAR -k KPOINTS and -i INCAR')
    parser.add_argument('-j', '--job', choices=['pchg','chg','md','nnff','mdnve','nnffnve', 'ini','zpe','mol','wav','opt','copt','sp','noD','kp','fake','neb','pseudo'], help='inquire for each file')
    parser.add_argument('-sj', '--subjob', default='sp', choices=['sp', 'cool', 'heat','quench'], help='sp for fake and others for md')
    parser.add_argument('-n', '--ndirs', default=5, type=int, help="number or dirs to make")
    ### POSCARs
    gposcar = parser.add_mutually_exclusive_group()
    gposcar.add_argument('-s', '--poscar', nargs='+', help='poscars in narrative mode')
    gposcar.add_argument('-p', '--prefix', help='select POSCAR using prefix lists for module common')
    ###
    parser.add_argument('-si', '--idposcar', nargs='+', help='in case poscar has index')
    parser.add_argument('-pot', '--potcar', choices=['new','potpaw-pbe-new','old','potpaw-pbe-old','potpaw-gga'], help='pseudo potential directory: ')
    parser.add_argument('-k', '--kpoints', nargs='+', help='input number of k-points in kx, ky, kz, or g for gamma')
    ### job dependent inputs
    #jobargs = parser.add_mutually_exclusive_group(title='job dependent arguments')
    parser.add_argument('-ja', '--jobadds', nargs='+', help='list of job dependent inputs')
    #jobargs.add_argument('-a', '--atoms', nargs='+', help='list of atoms for what?')
    #jobargs.add_argument('-sf', '--poscar_final',  help='for job=neb, input additional poscar')
    #jobargs.add_argument('-hpp', '--pseudoH', nargs='*', help='include pseudo H list ')
    ### KP tests
    g_ktest = parser.add_argument_group(title='KP tests')
    g_ktest.add_argument('-kd', '--kdim', default=3, type=int, choices=[1,2,3], help='input series of k-points [kx, ky, kz]*3')
    g_ktest.add_argument('-kps', '--kpoints_test', nargs='*', type=int, help='input series of k-points [kx, ky, kz]*3')
    ### toggle default: unset in the bare dir, set to j when INCAR.job exists
    parser.add_argument('-i', '--incar', help='in the order of designated_INCAR, INCAR.job,INCAR, dirname/INCAR]')
    parser.add_argument('-io', '--incar_option', nargs='*', help='to modify input INCAR]')
    parser.add_argument('-f', '--iofile', default='incar.key', help='only read file is possible')
    parser.add_argument('-d', '--dnames', nargs='+', help='get directory name')
    parser.add_argument('-al', '--all', action='store_true', help="prepared in job dir if not -s, -p, -k, -i")
    ### Running option
    g_run = parser.add_argument_group(title='Running options')
    g_run.add_argument('-ra', '--run_all', action='store_true', help="without asking")
    g_run.add_argument('-r', '--run', choices=['a','o','on','s','k'], help="o:overwrite run,on:overwrite stop,s:just submit,k:test input")

    ### VASP executable
    g_vasp  = parser.add_argument_group(title='VASP executable')
    g_vasp.add_argument('-exe', '--executable', choices=['gamma','xyrelax'], help='vasp execuatable: gamma, xy-relax')
    ### PBS arguments
    g_queue = parser.add_argument_group(title='QUEUE')
    g_queue.add_argument('-x', '--xpartition', type=int, help="partition in platinum")
    g_queue.add_argument('-N', '--nnode', type=int, help="number of nodes, can be used to calculate total nproc")
    g_queue.add_argument('-np', '--nproc', help="number of nproc, total for pt, per node for kisti ")
    g_queue.add_argument('-l', '--lkisti', nargs='*', help="kisti command line input")
    g_queue.add_argument('-o', '--option', help="error,exe; 'opt':converge, 'mem': lack, 'longnnn':long queue, 'g': gamma, 'vg':ML")
    args = parser.parse_args()

    ### running option
    if args.run_all:
        Lrun = 'a'
    elif args.run:
        Lrun = args.run
    else:
        Lrun = None

    pwd = os.getcwd()
    ### POSCARs and DIRECTORYs are abtained here and passed to make_vasp_dir()
    ### Apply dirnames to run fake job in KISTI
    job = args.job              # to pass job to function
    ### JOB = FAKE
    if args.job == 'fake':
        job = args.subjob
        ### not perfect
        if args.poscar:
            inposcar = args.poscar[0]
        else:
            if os.path.isfile('POSCAR'):
                inposcar = 'POSCAR'
            else:
                print("input poscar or make POSCAR in wdir")
                sys.exit(0)

        if args.dnames:
            if len(args.dnames) == 1:
                dirnames = []
                poscars  = []
                dname = args.dnames[0]
                ### make the poscars and args.dnames same for fake job in kisti 
                for a in list(string.ascii_lowercase)[:args.ndirs]:
                    dirnames.append(f"{dname}{a}")
                    # if
                    poscars.append(inposcar)
            else:
                poscars = args.dnames
        else:
            print("-j fake requires -d dnames")
            sys.exit(1)
    ### Normal JOBS: dirnames = list
    elif args.poscar:
        poscars=args.poscar
        if args.dnames:
            dirnames = args.dnames
        else:
            dirnames = get_dnames4pos(poscars)     # get dirname from POSCAR.name
    elif args.prefix:
        prefix0='POSCAR.'+args.prefix
        prefixes=[prefix0]
        print(f"{prefixes}")
        poscars = get_files_prefix(prefixes, pwd)
        if args.dnames:
            dirnames = args.dnames
        else:
            dirnames = get_dnames4pos(poscars)

    print(poscars)
    print(dirnames)
    if Lrun and 'k' in Lrun:
        print("stop before function")
        sys.exit(10)
    if len(poscars) != len(dirnames):
        print(f"poscar {len(poscars)} != directory {len(dirnames)}")
        sys.exit(11)
    #sys.exit(0)

    

    ### for kpoints-scan
    #if args.kp_test:
    if args.job == 'kp':
        if not args.lkisti:
            args.lkisti = 'kp'  # lkisti is changed to string from list
        kparray = np.array(args.kpoints_test)
        kps = kparray.reshape([-1,args.kdim])       # shape = n * 2 [ [1, 3], [2,2], [2, 3], ...] fir kdim=2
        for kp in kps:
            dname=[]
            if kp.size == 1:
                kp_in=[kp[0], kp[0], kp[0]]
            elif kp.size == 2:                      # [1, 3].size = 2
                kp_in=[kp[0], kp[1], 1]
            elif kp.size == 3:
                kp_in = list(kp)
            print(kp_in)
            kp_str = list(map(str, kp_in))
            dirname = dirnames[0]+'k'+list2str(kp_in)
            print(f"kp_string {kp_str}, dirname {dirname} in function {whereami()}()")
            dname.append(dirname)  # dname is string

            make_vasp_dir(job, poscars, args.potcar, args.jobadds, kp_str, args.incar, args.incar_option, dname, args.option, args.all, args.iofile, args.xpartition, args.nnode, args.nproc, args.executable, args.lkisti, Lrun)
    else:
################# 1       2        3       4                   5        6                7           8         9        10      11    12 13 14   15   16  
#def make_vasp_dir(job, subjob, poscars, apotcar, jobadds,           kpoints, incar,            dirnames, option, allprepared, iofile, qx,qN,qn,vasp_exe,lkisti,Lrun):
        make_vasp_dir(job, args.subjob, poscars, args.potcar, args.jobadds, args.kpoints, args.incar, args.incar_option, dirnames, args.option, args.all, args.iofile, args.xpartition, args.nnode, args.nproc, args.executable, args.lkisti, Lrun)
    return 0

if __name__ == '__main__':
    main()
