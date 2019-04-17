#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_files

def jobs(job,fname,que,np,nmem):
    if job == None:
        print("List this directory = ")
        mdir = os.path.dirname(__file__)
        exe, mod = dir_files(mdir)
        print("Executable:: ")
        exe_s = sorted(exe)
        mod_s = sorted(mod)
        for f in exe_s:
            print("    {}".format(f))
        print("Module:: ")
        for f in mod_s:
            print("    {}".format(f))
        print("#Comment for sge: -j [amp|grmx|qchem|sleep|qstat for q command|qsub for arguments]")
        print("    AMP::     sge_amp.tcsh")
        print("    Gromacs:: sge_mdrun2.tcsh")
        print("    Q-Chem::  sge_qchem.tcsh")
        print("    dummy job::  sge_sleep.csh")
    elif job == 'qstat':
        print("    $qstat -f    # see all nodes(/used/total) and my job")
        print("    $qfree       # free nodes")
        print("    $qhist       # shows number of running cores and jobs in Queue")
        print("    $qmem        # memory")
    elif job == 'qsub':
        print("    qsub -N jobname -v var=a_variable /qpfs/home/joonho/sandbox_gl/pypbs/sge_qchem.tcsh")
        print("\t$var can be used as variable in the script")
        print("\tscript name will be full name")
    elif job == 'qchem':
        print("    qsub -N opt -v job=1-PP-A-mul3-opt -v ver=\"3.2p\" /qpfs/home/joonho/sandbox_gl/pypbs/sge_qchem.tcsh")
        print("    in script::\n\t/gpfs/opt/qchem/bin/qchem -np 4 $job.in $job.out")
        print("\tenv variable such as $QC is not working in script")
    elif job == 'grmx':
        print("    qsub -v tpr=mdname sge_mdrun2.sh")
    elif job == 'amp':
        print("SGE validation::")
        print(f"    qsub -v fname={fname} -v py_job=val -v scan=scan -v nc={np} -q {que} -pe numa {np} -l h_vmem={nmem} $SGE_HOME/sge_amp.csh")
        print(      "\tpy_job for [validation||train||test]")
        print(      "\tscan variable (if exists) average for several trial with the same [HL, ene convergency]")
        print(      "\tsetenv SGE_HOME $HOME/sandbox_gl/pypbs")
        print(      "\tfind node:: sge_nodes.py || qstat -f")
        print("    in script::\n\tsetenv PYTHONPATH $HOME/sandbox_gl/pycommon:$HOME/sandbox_gl/mymplot")
        print(                 "\tset PYTHON = \"$HOME/anaconda3/bin/python\"")
        print(                 "\tset EXE = \"$HOME/sandbox_gl/py_ai/amp_ene.py\"")
        print(                 "\t$PYTHON $EXE $fname $py_job -hl 4 4 4 -el 0.001 -n 5 -g")
    elif job == 'qsleep':
        print("Dummy job::")
        print(f"            qsub -q {que} -pe numa {np} -l mem={nmem} $SB/pypbs/sge_sleep.csh")
        print("    Sample: qsub -q qname(sandy@slet02) -pe numa 6(num of process) -l mem=24G $SB/pypbs/sge_sleep.csh")

    return
        

def main():

    parser = argparse.ArgumentParser(description="display Usage for /pycommon in SGE  ")
    parser.add_argument('-j','--job',  choices=['amp','grmx','qchem','qstat','qsub','qsleep'], help="several explanation option ")
    #parser.add_argument('-l','--list', action='store_true',  help="list directory files ")
    parser.add_argument('-f', '--file', default=".extxyz", help='energy data input file')
    parser.add_argument('-q', '--queue', default="sandy@opt03", help='if you want to assign queue')
    parser.add_argument('-n', '--nprocess', default=12, type=int,  help='number of parallel process')
    parser.add_argument('-l', '--memory', default=2.0, type=float,  help='number of parallel process')
    args = parser.parse_args()

    jobs(args.job, args.file, args.queue, args.nprocess, args.memory)

if __name__ == "__main__":
    main()
