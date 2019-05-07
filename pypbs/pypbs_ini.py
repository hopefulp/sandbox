#!/gpfs/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_files, yes_or_no

def jobs(work,fname,Lrun,quejob,np,nmem,que):
    if work == None:
        mdir = os.path.dirname(__file__)
        print(f"List directory in {mdir}")
        exe, mod = dir_files(mdir)
        print("Executable:: ")
        exe_s = sorted(exe)
        mod_s = sorted(mod)
        for f in exe_s:
            print("    {}".format(f))
        print("Module:: ")
        for f in mod_s:
            print("    {}".format(f))
        print("#Comment for sge: -j [sge|amp|grmx|qchem|sleep|qsub for arguments]")
        print("    AMP::     sge_amp.tcsh")
        print("    Gromacs:: sge_mdrun2.tcsh")
        print("    Q-Chem::  sge_qchem.csh")
        print("    dummy job::  sge_sleep.csh")
    elif work == 'sge':
        print("    $qstat -f    # see all nodes(/used/total) and my job")
        print("    $qstatf      # aliast qstat -f ...")
        print("    $qfree       # free nodes")
        print("    $qhist       # shows number of running cores and jobs in Queue")
        print("    $qmem        # memory")
    elif work == 'qsub':
        print("    qsub -N jobname -v var=a_variable /gpfs/home/joonho/sandbox_gl/pypbs/sge_qchem.csh")
        print("\t$var can be used as variable in the script")
        print("\tscript name will be full name")
    elif work == 'qchem':
        ### version control
        #com = f"    qsub -N {quejob} -v qcjob={fname} -v ver=\"5.1p\" -v nc={np} -pe numa {np} /gpfs/home/joonho/sandbox_gl/pypbs/sge_qchem.csh"
        if fname == None:
            com = f"    qsub -N {quejob} -v qcjob=a.in -v nc={np} -pe numa {np} /gpfs/home/joonho/sandbox_gl/pypbs/sge_qchem.csh"
        else:
            com = f"    qsub -N {quejob} -v qcjob={fname} -v nc={np} -pe numa {np} /gpfs/home/joonho/sandbox_gl/pypbs/sge_qchem.csh"
        print(com)
        if not Lrun:
            print("\tuse -r to queue submit")
            print("    \"{}\" is aliased to \"{}\"".format(os.path.basename(__file__)+" -w qchem", "qchem_run"))
            print("    use :: qchem_run -j Queue-jobname -f qin-file -n num_process")
            print("    e.g.:: qchem_run -jopt -fCO2 -n4")
            print("For batch::\n\tdir_run.sh in \# $1 for qchem input type")
            if fname == None:
                print(f"    in script::\n\t/gpfs/home/joonho/sciware/qchem5.1p/bin/qchem -np {np} a.in a.out")
            else:
                print(f"    in script::\n\t/gpfs/home/joonho/sciware/qchem5.1p/bin/qchem -np {np} {fname}.in {fname}.out")
            print("\tenv variable such as $QC is not working in script")
        else:
            q = "will you run?"
            if yes_or_no(q):
                os.system(com)
    elif work == 'grmx':
        print("    qsub -v tpr=mdname sge_mdrun2.sh")
    elif work == 'amp':
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
    elif work == 'qsleep':
        print("Dummy job::")
        print(f"            qsub -q {que} -pe numa {np} -l mem={nmem} $SB/pypbs/sge_sleep.csh")
        print("    Sample: qsub -q qname(sandy@slet02) -pe numa 6(num of process) -l mem=24G $SB/pypbs/sge_sleep.csh")

    return
        

def main():

    parser = argparse.ArgumentParser(description="display Usage for /pycommon in SGE  ")
    parser.add_argument('-w','--work',  choices=['amp','qchem','sge','qsub','qsleep','grmx'], help="several explanation option ")
    parser.add_argument('-f', '--file', help='energy data input file')
    parser.add_argument('-r', '--run', action='store_true', help='run command')
    queue = parser.add_argument_group(title = 'SGE Queue')
    queue.add_argument('-j', '--que_job', default="opt", help='SGE Queue jobname')
    queue.add_argument('-q', '--queue', default="sandy@opt03", help='if you want to assign queue')
    queue.add_argument('-n', '--nprocess', default=12, type=int,  help='number of parallel process')
    queue.add_argument('-l', '--memory', default=2.0, type=float,  help='size of memory per process')
    ### using add_argument_group, if develop other works, try add_subparsers
    #qchem = parser.add_argument_group(title='Q-Chem')
    args = parser.parse_args()

    jobs(args.work,args.file,args.run,args.que_job,args.nprocess,args.memory,args.queue)

if __name__ == "__main__":
    main()
