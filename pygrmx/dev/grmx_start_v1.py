#!/usr/bin/python
import argparse
import os
import common
from shutil import copy

def grom_pp(fpre, mfile, mdyn,fgro):
    print("Do Grmx preprocess in a given system")
    ftarget = fpre + "_" + mdyn
    cmd = "grompp -f %s -c %s -p %s.top -o %s.tpr" % (mfile,fgro,fpre,ftarget)
    print(cmd)
    print("-f input .mdp")
    print("-c input .gro, grmx coordinate")
    print("-p input .top")
    print("-o output .tpr for MDrun")
    q2 = "will you run ? "
    if common.yes_or_no(q2):
        os.system(cmd)      # this does not make .tpr
        if os.path.isfile('%s.tpr' % ftarget):
            print("%s.tpr was made" % ftarget)
            return 1
        else:
            print("Fail to generate %s.tpr" % ftarget)
    else:
        print("Job terminated")
    return 0

def run_job(job, ifile, icharge, mdfile, dyn, jobname, grofile):
    pdir = os.getcwd()
    if not os.path.isfile(ifile):
        print("there is no %s" % ifile)
        exit(10)
    f_pre, f_ext = common.fname_decom(ifile)
    L_write=False
    q = "continue to overwrite ? "
    q2 = "will you run ? "
    ### step 1: make .pdb file
    if job == 'mol2pdb':
        if os.path.isfile('%s.pdb' % f_pre):
            print("%s.pdb exists: " % f_pre)
            if common.yes_or_no(q):
                L_write=True
        else:
            L_write=True
        if L_write:
            com = "mol2pdb.pl %s.in > %s.pdb " % (f_pre, f_pre)
            print(com)
            q1 = "will you make .pdb from .in(q-chem) ?"
            print("input: \"a.in\"\noutput: \"a.pdb\"")
            if common.yes_or_no(q1):
                os.system(com)
                print("Step 1: %s.pdb was made" % f_pre)
                print("continue to Step 2: pdb2mol2")
            else:
                print("job was terminated")
            return 0
    ### step 2: make .mol2 file
    elif job == 'pdb2mol2':
        if os.path.isfile('%s.mol2' % f_pre):
            print("%s.mol2 exists: " % f_pre)
            if common.yes_or_no(q):
                L_write=True
        else:
            L_write=True
        if L_write:
            print("input: \"a.pdb\"\noutput: \"a.mol2\"")
            print("net charge: %d" % icharge)
            print("atom type: gaff")
            print("to modify more, correct source code")
            com = "antechamber -i %s.pdb -fi pdb -o %s.mol2 -fo mol2 -at gaff -c bcc -nc %d" % (f_pre, f_pre, icharge)
            print("%s \n This job will take time" % com)
            q1 = "will you make .mol2 from .pdb ?"
            if common.yes_or_no(q1):
                os.system(com)
                print("%s.mol2 was made" % f_pre)
                print("continue to make frcmod step 3: frcmod")
            else:
                print("job was terminated")
            return 0
    ### Step 3: make .frcmod file
    elif job == 'frcmod':
        if os.path.isfile('%s.frcmod' % f_pre):
            print("%s.frcmod exists: " % f_pre)
            if common.yes_or_no(q):
                L_write=True
        else:
            L_write=True
        if L_write:            
            com = "parmchk2 -i %s.mol2 -f mol2 -s 1 -o %s.frcmod" % (f_pre, f_pre)
            print(com)
            print("-i for input file")
            print("-f for input format")
            print("-s for FF parameter set")
            print("-o for frcmod output file")
            q1 = "will you make .frcmod from .mol2 ?"
            if common.yes_or_no(q1):
                os.system(com)
                print("%s.frcmod was made" % f_pre)
                print("continue to check mol2 charge w mol2_chg")
            else:
                print("job was terminated")
            return 0
    ### Step 4: check .mol2-charge then go to tleap
    elif job == "mol2chg":
        com = 'cal_chg.sh %s' % ifile
        print("check total charge of %s" % ifile)
        os.system(com)
        q1 = "is charge OK ?"
        if common.yes_or_no(q1):
            print("Alert:: tleap is interactive mode, Do this on shell")
            print("$tleap")
            print("tleap> source leaprc.gaff")
            print("tleap> mod = loadamberparams %s.frcmod" % f_pre)
            print("tleap> x = loadmol2 %s.mol2" % f_pre)
            print("tleap> saveamberparm x %s.top %s.crd" % (f_pre, f_pre))
            print("tleap> quit")
            print("with %s.top %s.crd, go to step 5: amber2grmx" % (f_pre, f_pre))
        else:
            print("job was terminated")
            print("modify %s.mol2 then go to Step 4: mol2_chg again")
        return 0
    ### step 5: tleap-made .top and .crd. make gmx.top, gmx.gro        
    elif job == "amber2grmx":
        q1 = "top, crd were made from tleap ? "
        if common.yes_or_no(q1):
            if os.path.isfile("%s.top" % f_pre) and os.path.isfile("%s.crd" % f_pre):
                com = "acpype.py -p %s.top -x %s.crd" % (f_pre, f_pre)
                print(com)
                print("-i for input but not used here")
                print("-p amber parameter topology, always used with -x")
                print("-x amber input coordinate, always used with -p")
                os.system(com)
                print("_GMX.top, _GMX.gro were made")
            else:
                print("this work is not ready")
                exit(44)
        q1 = "will you arrange directory ? "
        if common.yes_or_no(q1):
            tmpdir = pdir + "/qc_input"
            if not os.path.isdir(tmpdir):
                os.mkdir(tmpdir)
            os.system("mv %s/* %s" % (pdir, tmpdir))
            os.chdir(tmpdir)
            os.system("cp *GMX* %s" % pdir)
            os.chdir(pdir)
            print("change gro and top file to a proper names")
            print("modify .gro for box size: use ins/substitute mode to keep format")
            print("editconf to move the system to center if necessary")
            print("modify .top charge again: got to Step 5.1 \"topchg\"")
        return 0            
    ### Step 5.1: check gromax.top charge then continue
    elif job == "topchg":
        print("check gromac %s.top charge" % f_pre)
        com = "cal_chg.sh %s" % ifile
        os.system(com)
        print("check or modify %s.top file then go to Step 6 to make water box: waterbox" % f_pre)
        return 0
    ### step 6: make water box in gromacs .top .gro            
    elif job == "waterbox":
        print("it will fill water in the box")
        com = "genbox -cp %s.gro -cs spc216.gro -o %s_sol.gro -p %s.top" %(f_pre, f_pre, f_pre)
        print(com)
        print("-cp input for gro structure file with proper box size")
        print("-cs input for sovent water structure")
        print("-o output for _sol.gro")
        print("-p in/out topology file is overwritten")
        q1 = "will you run ? "
        if common.yes_or_no(q1):
            os.system(com)
            print("%s_sol.gro was made" % f_pre)
            print("%s.top was overwritten" % f_pre)
            print("modify .top for water, Step 7 : top_gaff_addwater")
        else:
            print("The job was terminated")
        return 0
    ### Step 7: modify .top for water before making pp
    elif job == "modtop_water":
        print("modify .top file on the fly")
        if f_ext != "top":
            print("input .top")
            exit(55)
        cmd = "top_gaff_add_water.sh %s.top" % f_pre
        os.system(cmd)
        print("Continue to step 8 with %s_sol.gro: pre_pp" % f_pre)
        return 0
    ### step 8: make MD inputfile of .tpr for system without anion
    elif job == "pre_pp":
        if f_ext != "top":
            print("input .top")
            exit(56)
        if not mdfile:
            print("input mdfile using --md")
            exit(57)
        print("Do Grmx preprocess in a given system")
        if grofile:
            gfile = grofile
        else:
            gfile = f_pre + "_sol.gro"
        Lsuccess = grom_pp(f_pre, mdfile, dyn, gfile)
        if Lsuccess:
            print("Continue to modify .top before adding anion to system, step 9: top_gaff_addion")
        return 0
    ### stop 9: add anion to MD input file
    elif job == "modtop_ion":
        if f_ext != "top":
            print("input .top")
            exit(55)
        print("modify .top file on the fly for adding ion")
        cmd = "top_gaff_add_ion.sh %s.top" % f_pre
        print(cmd)
        os.system(cmd)
        print("%s.top was modified" % f_pre)
        print("-- Add ion interactively: -- ")
        cmd = "genion -s %s.tpr -o %s_sol_neut.gro -p %s.top -nn 1 -nname CL -nq -1" % (f_pre,f_pre, f_pre)
        print(cmd)
        print("-s input of MD running file")
        print("-o output .gro of coordination")
        print("-p in/out .top overwritten")
        q1 = "will you run interactive commander genion ?"
        if common.yes_or_no(q1):
            os.system(cmd)
        print("after made %s_sol_neut.gro, make .tpr again: go to step 10 using make_pp" % f_pre)
        return 0 
    ### step 10: make pp
    if job ==  'do_pp':
        if f_ext != "top":
            print("input .top")
            exit(96)
        if not mdfile:
            print("input mdfile using --md")
            exit(97)
        if grofile:
           gfile = grofile
        else:
            gfile = f_pre + "_sol_neut.gro"
        Lsuccess = grom_pp(f_pre, mdfile, dyn, gfile)
        if Lsuccess:
            print("Complete %s_%s.tpr was made: ready to run MD" % (f_pre, dyn) )
            print("submit a job with -v variable, 'submit': qsub -N jobname -v tpr=jobname sge_mdrun.tcsh")
            print("note: use .tcsh: bash script .sh cannot be used")
            print("Or submit at last step:: 'submit'")
        return 0
    ### Step last: submit a job 
    if job == 'submit':
        home = os.environ['HOME']
        pbs = home + '/sandbox_gl/pypbs/sge_mdrun.tcsh'
        #print(home)
        ### jobname makes directory
        tgro = pdir + '/' + f_pre + ".gro"
        if f_ext != "tpr":
            print("run with .tpr name")
            exit(99)
        elif os.path.isfile(tgro):
            print("do not overwrite %s.gro" % f_pre)
            exit(98)
            
        if jobname:
            jname = jobname
            job_dir = pdir + '/' + jname
            os.mkdir(job_dir)
            fname = pdir + '/' + f_pre + '.tpr'
            copy(fname, job_dir)
            os.chdir(job_dir)
        else:
            #jname = f_pre
            print("try job name with -n or --jobname")
            exit(100)
        cmd = "qsub -N %s -v tpr=%s %s" % (jname, f_pre, pbs)
        print(cmd)
        print("job: mdrun -s %s.tpr -c %s.gro -o %s.trr -e %s.edr -g %s.log" %(f_pre, f_pre, f_pre, f_pre, f_pre))
        print("input: %s.tpr" % f_pre)
        print("output: %s.gro, %s.trr, %s.edr, %s.log" % (f_pre, f_pre, f_pre, f_pre))
        q = 'submit this job ?'
        if common.yes_or_no(q):
            os.system(cmd)
            print("job submited.")
        else:
            print("job skipped")
        return 0

def main():
    parser = argparse.ArgumentParser(description='run MD-Gromacs step by step w. GAFF; for last step of grompp and mdrun, go to step 10 of "do_pp"')
    parser.add_argument('work', choices=['mol2pdb','pdb2mol2','frcmod','mol2chg','amber2grmx','topchg','waterbox','modtop_water','pre_pp','modtop_ion','do_pp','submit'], help='individual jobs from qchem input file to gromacs MD')
    parser.add_argument('-i', '--infile', help='always input file is expected')
    parser.add_argument('-c', '--charge', default=0, type=int, help="total charge is for pdb2mol2")
    parser.add_argument('-f', '--mdfile', help="MD parameter file")
    parser.add_argument('-d', '--dynamics', default='md', choices=['md', 'min'], help="MD parameter file")
    parser.add_argument('-j', '--jobname', help="input job name in qsub")
    parser.add_argument('-u', '--usage', action='store_true', help=( "Usage of each work :: %s  work-type -u" % os.path.basename(__file__)))
    parser.add_argument('-g', '--grofile', help="gro file used for grompp")
    args = parser.parse_args()

    if args.usage:
        #print("Usage:: ", end='')
        #print("work list = mol2pdb, pdb2mol2, frcmod, mol2chg, amber2grmx, topchg, waterbox, top_gaff_addwater, pre_pp, top_gaff_addion, do_pp, submit")
        print "Usage :: ",
        if args.work == 'mol2pdb':
            print("%s %s -i a.in" % (os.path.basename(__file__), args.work))
            print("\tmakes a.pdb")
        elif args.work == 'pdb2mol2':
            print("%s %s -i a.pdb -c \"int\" " % (os.path.basename(__file__), args.work))
            print("\tmakes a.mol2")
        elif args.work == 'frcmod':
            print("%s $s -i a.mol2" % (os.path.basename(__file__), args.work))
            print("\tmakes a.frcmod")
        elif args.work == 'mol2chg':
            print("%s %s -i a.mol2" %  (os.path.basename(__file__), args.work))
            print("Change charge of a.mol2 manually before advance to tleap")
        elif args.work == 'amber2grmx':
            print("%s %s -i a.top" %  (os.path.basename(__file__), args.work))
            print("\tmakes gromacs .top and .gro")
        elif args.work == 'topchg':
            print("%s %s -i a.top" %  (os.path.basename(__file__), args.work))
            print("\tcheck .top charge before water box: correct .top chg and .gro box-size")
        elif args.work == 'waterbox':
            print("%s %s -i a.top" %  (os.path.basename(__file__), args.work))
            print("\tfill box with water")
        elif args.work == 'modtop_water':
            print("%s %s -i a.top" %  (os.path.basename(__file__), args.work))
            print("\tmodify a.top for water FF file")
        elif args.work == 'pre_pp':
            print("%s %s -i a.top -f mdpfiile.mdp -d [md|dy]" %  (os.path.basename(__file__), args.work))
            print("\tmakes .tpr MD running file")
        elif args.work == "modtop_ion":
            print("%s %s -i a.top" %  (os.path.basename(__file__), args.work))
            print("\tmakes new .gro and overwrite .top")
        elif args.work == "do_pp":
            print("%s %s -i a.top -f mdpfiile.mdp -d [md|dy] [-g grofile]" %  (os.path.basename(__file__), args.work))
            print("\tmakes the complete .tpr MD running file: -d for fname.tpr")
        elif args.work == "submit":
            print("%s %s -i a.tpr -j jobname" % (os.path.basename(__file__), args.work))
            print("\tmake dir[jobname] for dynamics and submit job there")
        return 0

    #run_job(args.work, args.infile, args.charge, args.mdfile, args.dynamics, args.jobname)
    run_job(args.work, args.infile, args.charge, args.mdfile, args.dynamics, args.jobname,args.grofile)
    return 0

if __name__ == '__main__':
    main()
