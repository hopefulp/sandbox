#!/gpfs/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_files

def jobs(job):
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
        print("#Comment for sge: -j [amp|grmx|qchem|qstat for q command|qsub for arguments]")
        print("    sge_amp.tcsh for AMP")
        print("    sge_mdrun2.tcsh for GROMACS MD")
        print("    sge_qchem.tcsh for Q-Chem")
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
        print("    qsub -v fname=water128.extxyz -v py_job=tr /gpfs/home/joonho/sandbox_gl/pypbs/sge_amp.tcsh")
        print("    in script::\n\t/gpfs/home/joonho/sandbox_gl/py_ai/sge_amp_ene.py $fname $py_job -hl 4 4 4 -el 0.001 +g")

    return
        

def main():

    parser = argparse.ArgumentParser(description="display Usage for /pycommon in SGE  ")
    parser.add_argument('-j','--job',  choices=['amp','grmx','qchem','qstat','qsub'], help="several explanation option ")
    #parser.add_argument('-l','--list', action='store_true',  help="list directory files ")
    args = parser.parse_args()

    jobs(args.job)

if __name__ == "__main__":
    main()
