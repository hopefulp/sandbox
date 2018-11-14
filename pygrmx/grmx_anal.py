#!/usr/bin/python
import argparse
import os
import common
from shutil import copy

def run_job(job, ifile, eoption, nskip, pbctype):
    pdir = os.getcwd()
    if not os.path.isfile(ifile):
        print("there is no %s" % ifile)
        exit(10)
    f_pre, f_ext = common.fname_decom(ifile)
    L_write=False
    q = "overwrite ? "
    if job == 'ene':
        if os.path.isfile('%s.xvg' % f_pre):
            print("%s.xvg exists: " % f_pre)
            if common.yes_or_no(q):
                L_write=True
        else:
            L_write=True
        if L_write:
            com = "g_energy -f %s.edr -dp -o %s.xvg " % (f_pre, f_pre)
            print(com)
            os.system(com)
            if os.path.isfile('%s.xvg' % f_pre):
                print("%s.xvg was made" % f_pre)
            else:
                print("job failed")
        return 0            
    elif job == 'trj':
        if os.path.isfile('%s.gro' % f_pre):
            print("%s.gro exists: " % f_pre)
            if common.yes_or_no(q):
                L_write=True
        else:
            L_write=True
        if L_write:
            com = "trjconv -f %s.trr -s %s.tpr -o %s.gro -skip %d -pbc %s"  % (f_pre, f_pre, f_pre, nskip,pbctype)
            print(com)
            q2="will you run?"
            if common.yes_or_no(q2):
                os.system(com)
                print("%s.gro was made" % f_pre)
            else:
                print("job aborted")
        return 0

def main():
    parser = argparse.ArgumentParser(description='run MD-Gromacs step by step w. GAFF')
    parser.add_argument('job', choices=['ene', 'trj'], help='analysis of MD results')
    parser.add_argument('inp_file', help='always input file is expected')
    parser.add_argument('-e', '--eoption', default='tE', help='input diverse energy option')
    parser.add_argument('-s', '--skip', type=int, default=10000, help='get trajectories in every skip steps')
    parser.add_argument('-p', '--pbc', default="mol", help='for molecular connection in vmd')

    args = parser.parse_args()

    run_job(args.job, args.inp_file, args.eoption, args.skip, args.pbc)
    return 0

if __name__ == '__main__':
    main()
