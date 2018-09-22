#!/usr/bin/python
import argparse
import os
import common

def run_job(job, ifile, icharge, mdfile):
    if not os.path.isfile(ifile):
        print("there is no %s" % ifile)
        exit(10)
    
    f_pre, f_ext = common.fname_decom(ifile)
    md_pre, md_ext = common.fname_decom(mdfile)
    L_write=False
    q = "overwrite ? "
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
            os.system(com)
            print("%s.pdb was made" % f_pre)
            return 0

    elif job == 'pdb2mol2':
        if os.path.isfile('%s.mol2' % f_pre):
            print("%s.mol2 exists: " % f_pre)
            if common.yes_or_no(q):
                L_write=True
        else:
            L_write=True
        if L_write:            
            com = "antechamber -i %s.pdb -fi pdb -o %s.mol2 -fo mol2 -at gaff -c bcc -nc %d" % (f_pre, f_pre, icharge)
            print("%s will take time" % com)
            print("-at : atom type")
            print("-c  : charge calculation method")
            print("-nc : net charge of the system")
            #os.system(com)
            print("%s.mol2 was made" % f_pre)
            q2 = "make frcmod file ?"
            if common.yes_or_no(q2):
                com = "parmchk2 -i %s.mol2 -f mol2 -s 1 -o %s.frcmod" % (f_pre, f_pre)
                print(com)
                os.system(com)
                print("%s.frcmod was made" % f_pre)
            return 0
    elif job == "mol2_chg":
        com = 'mol2_chg.sh %s' % ifile
        print("check total charge of %s" % ifile)
        os.system(com)
        print("if ok go to 'tleap'")
        return 0
    elif job == "amber2grmx":
        q = "top, crd were made from tleap ? "
        if common.yes_or_no(q):
            if os.path.isfile("%s.top" % f_pre) and os.path.isfile("%s.crd" % f_pre):
                com = "acpype.py -p %s.top -x %s.crd" % (f_pre, f_pre)
                print(com)
                os.system(com)
                print("_GMX.top, _GMX.gro were made")
        q = "will you arrange directory ? "
        if common.yes_or_no(q):
            pdir = os.getcwd()
            wdir = pdir + "/tmp"
            if not os.path.isdir(wdir):
                os.mkdir(wdir)
            os.system("mv %s/* %s" % (pdir, wdir))
            os.chdir(wdir)
            os.system("cp *GMX* %s" % pdir)
            os.chdir(pdir)
            print("change gro and top file to a proper names")
            print("modify .gro for box size")
            print("editconf to move the system to center if necessary")
            print("modify .top for syntax")
    elif job == "waterbox":
        print("it will fill water in the box")
        com = "genbox -cp %s.gro -cs spc216.gro -o %s_sol.gro -p %s.top" %(f_pre, f_pre, f_pre)
        print(com)
        #os.system(com)
        print("%s_sol.gro was made" % f_pre)
        print("%s.top was overwritten" % f_pre)
        return 0
    elif job == "top_gaff_addwater":
        print("modify .top file on the fly")
        if f_ext != "top":
            print("input .top")
            exit(55)
        cmd = "top_gaff_add_water.sh %s.top" % f_pre
        os.system(com)
        return 0
    elif job == "pre_pp":
        if f_ext != "top":
            print("input .top")
            exit(56) 
        print("Do Grmx preprocess in a given system")
        cmd = "grompp -f %s -c %s_sol.gro -p %s.top -o %s.tpr" % (mdfile,f_pre,f_pre,f_pre)
        print(cmd)
        os.system(cmd)
        print("%s.tpr was made" % f_pre)
        return 0
    elif job == "top_gaff_addion":
        if f_ext != "top":
            print("input .top")
            exit(55)
        print("modify .top file on the fly for adding ion")
        cmd = "top_gaff_add_ion.sh %s.top" % f_pre
        print(cmd)
        os.system(cmd)
        print("%s.top was modified" % f_pre)
        print("for %s.top was modified, add ion interactively in %s_sol_neut.gro and %s.top" % (f_pre, f_pre, f_pre))
        cmd = "genion -s %s.tpr -o %s_sol_neut.gro -p %s.top -nn 1 -nname CL -nq -1" % (f_pre,f_pre, f_pre)
        print("use :: %s" % cmd)
        print("%s_sol_neut.gro will be made and %s.top will be overwritten" % (f_pre,f_pre))
        return 0 
    elif job == "gen_pp":
        print("if complete %s.gro, %s.top was made " % (f_pre, f_pre))
        cmd = "grompp -f %s -c %s_sol_neut.gro -p %s.top -o %s_%s.tpr" % (mdfile, f_pre, f_pre, f_pre, md_pre)
        print(cmd)
        os.system(cmd)
        print("%s_%s.tpr was newly generated" % (f_pre, md_pre))
        print("run mdrun")
        return 0
    elif job == "md_run":
        pass
        return 0


def main():
    parser = argparse.ArgumentParser(description='run MD-Gromacs step by step w. GAFF')
    parser.add_argument('job', choices=['mol2pdb', 'pdb2mol2', 'mol2_chg', 'amber2grmx', 'waterbox', 'top_gaff_addwater', 'pre_pp','top_gaff_addion', 'gen_pp', 'md_run'], help='individual jobs from qchem input file to gromacs MD')
    parser.add_argument('inp_file', help='always input file is expected')
    parser.add_argument('-c', '--charge', default=0, type=int, help="total charge is for pdb2mol2")
    parser.add_argument('-d', '--md', default='md.mdp', help="MD parameter file")
    args = parser.parse_args()

    run_job(args.job, args.inp_file, args.charge, args.md)
    return 0

if __name__ == '__main__':
    main()
