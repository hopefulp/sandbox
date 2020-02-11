#!/home/joonho/anaconda3/bin/python

import argparse
import os
import sys
from common import yes_or_no, list2str

amp_db = ['amp-fingerprint-primes.ampdb', 'amp-neighborlists.ampdb', 'amp-fingerprints.ampdb']

def make_dir(HL, elimit):
    if HL:
        dname = "".join(str(HL))        # list to string
        print(dname)
        dname += str(elimit)
        return "HL"+dname
    else:
        print("HL is empty")
        sys.exit(10)

def submit(qjobname,ncore,infile,amp_job,nHL,elimit):
    pwd = os.getcwd()
    str_hl = " ".join(str(x) for x in nHL)
    #if yes_or_no("run ?"):
    #jobdir = make_dir(nHL, elimit)
    jobdir = 'HL'+list2str(nHL)+str(elimit)

    if jobdir in os.listdir(cwd):
        print(f"{jobdir} exists")
        return 0
    else:
        com = f"will you make {jobdir}?"
        if yes_or_no(com):
            os.mkdir(jobdir)
            print(f"{jobdir} was made")
            new_dir = pwd + "/" + jobdir
            #os.system(f"cp -r amp-* {new_dir}")
            os.chdir(f"{new_dir}")
            for amp_dir in amp_db:
                os.system(f"ln -s {pwd}/{amp_dir} {new_dir}/{amp_dir}")
            os.system(f"cp {pwd}/OUTCAR {new_dir}")
            print("amp- directories and OUTCAR were copied")
        comm = f"qsub -N {qjobname} -pe numa {ncore} -v fname={infile} -v pyjob={amp_job} -v nc={ncore} -v hl=\"{str_hl}\" -v el={elimit} $SB/pypbs/sge_amp.csh"
        str1 = f"will you run: \n{comm}"
        if yes_or_no(str1):
            os.system(comm)
        return 0

def main():
    parser = argparse.ArgumentParser(description="submit jobs in queue in SGE(mlet) ")
    parser.add_argument('-qj', '--queue_jobname', default='amp', help='job option:"train","test","md","validation","profile"')
    parser.add_argument('-f', '--input_file', default='OUTCAR', help='ASE readible file: extxyz, OUTCAR(VASP) ')
    parser.add_argument('-j', '--amp_job', default='tr', help='job option:"train","test","md","validation","profile"')
    parser.add_argument('-hl', '--hidden_layers', nargs='*', type=int, default=[8,8,8], help='Hidden Layer of lists of integer')
    parser.add_argument('-nc', '--ncore', default=16, type=int, help='number of core needs to be defined')
    parser.add_argument('-el', '--e_convergence', default=0.0001, type=float, help='energy convergence limit')
    args = parser.parse_args()
    print(args.hidden_layers)
    submit(args.queue_jobname, args.ncore, args.input_file, args.amp_job, args.hidden_layers, args.e_convergence)

if __name__ == "__main__":
    main()
