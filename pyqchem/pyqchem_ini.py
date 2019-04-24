#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_files

usage = {   'xyz22mol' : ' a.xyz[a.mol]\t# makes a.mol[a.xyz]',
            'qout_geo' : ' a.out\t# extract the last optimized geometry as a.xyz',
            'qout_non_opt': ' a.out\t# extract the last non-optimized geometry as a.xyz'
        }
def if_usage(f):
    f_pre = f.split('.')[0]
    if f_pre in usage.keys():
        print("    {}".format(f), usage[f_pre])
    else:
        print("    {}".format(f))
def jobs(job,ifile,np):
    if job == None or re.search("cl", job):
        print("List this directory = ")
        mdir = os.path.dirname(__file__)
        exe, mod = dir_files(mdir)
        print("Executable:: ")
        sort_exe = sorted(exe)
        sort_mod = sorted(mod)
        if job == None:
            for f in sort_exe:
                print("    {}".format(f))
        else:
            lxyz=[]
            qc_out=[]
            for f in sort_exe:
                if re.match('xyz',f):
                    lxyz.append(f)
                    continue
                elif re.match('qout',f):
                    qc_out.append(f)
                    continue
                print("    {}".format(f))
            ### classify xyz files
            print("  {:<10}::".format('XYZ format'))
            for f in lxyz:
                if_usage(f)
            ### classify qout files
            print("  {:<10}::".format('QC-out'))
            for f in qc_out:
                if_usage(f)
        print("Module:: ")
        for f in sort_mod:
            print("    {}".format(f))
        print("#Comment: -j run for 'how to run'\n\t -j classify for detail ")
    elif job == 'run':
        com_serial="qchem {0}.in {0}.out &".format(ifile)
        com_parallel="mpirun -np {1} $QC/exe/qcprog {0}.in $QCSCRATCH > {0}.out &".format(ifile,np)
        print("{:^8}::".format('Chi'))
        print("\t{:^8}::".format('serial'), com_serial)
        print("\t{:^8}::".format('parallel'), com_parallel)

    

def main():

    parser = argparse.ArgumentParser(description="explanation for /pyqchem ")
    parser.add_argument('-j','--job', choices=['run','classify'],  help="qchem run in chi, mlet ")
    parser.add_argument('-f','--infile',  help="qchem input file")
    parser.add_argument('-np','--nprocess', default=2, type=int, help="number of parallel process")
    args = parser.parse_args()

    jobs(args.job,args.infile, args.nprocess)

if __name__ == "__main__":
    main()
