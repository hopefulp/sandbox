#!/gpfs/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_all, MyClass, dir_classify

sge     = MyClass()
grmx    = MyClass()
qsub    = MyClass()
usage   = MyClass()
amp     = MyClass()
qchem   = MyClass()
qsleep  = MyClass()

# if not use input variables, define here
sge.usage="$qstat -f    # see all nodes(/used/total) and my job\
            \n\t\t$qstatf      # aliast qstat -f ...\
            \n\t\t$qfree       # free nodes\
            \n\t\t$qfree       # free nodes\
            \n\t\t$qhist       # shows number of running cores and jobs in Queue\
            \n\t\t$qmem        # memory\
            "
qsub.usage="qsub -N jobname -v var=a_variable /gpfs/home/joonho/sandbox_gl/pypbs/sge_qchem.csh\
            \n\t\t$var can be used as variable in the script\
            \n\t\tscript name should be full name\
            "
grmx.usage="qsub -v tpr=mdname sge_mdrun2.sh"

usage.gen="#Comment for sge: -w [sge|amp|grmx|qchem|sleep|qsub for arguments]\
            \n\t\tAMP::     sge_amp.tcsh\
            \n\t\tGromacs:: sge_mdrun2.tcsh\
            \n\t\tQ-Chem::  sge_qchem.csh\
            \n\t\tdummy job::  sge_sleep.csh"
amp.usage="\n\t\tpy_job for [validation||train||test]\
            \n\t\tscan variable (if exists) average for several trial with the same [HL, ene convergency]\
            \n\t\tsetenv SGE_HOME $HOME/sandbox_gl/pypbs\
            \n\t\tfind node:: sge_nodes.py || qstat -f\
            \n\tin script::\n\tsetenv PYTHONPATH $HOME/sandbox_gl/pycommon:$HOME/sandbox_gl/mymplot\
            \n\t\tset PYTHON = \"$HOME/anaconda3/bin/python\"\
            \n\t\tset EXE = \"$HOME/sandbox_gl/py_ai/amp_ene.py\"\
            \n\t\t$PYTHON $EXE $fname $py_job -hl 4 4 4 -el 0.001 -n 5 -g\
            "
qchem.usage=""
qsleep.usage=""
classobj_dict={'SGE': sge, 'GRMX': grmx, 'AMP':amp, 'USAGE': usage, 'Qsub': qsub, 'QChem': qchem, 'Qsleep': qsleep}
lclass_var=['QChem', 'AMP', 'Qsleep']

def works(work,Lusage,fname,Lrun,quejob,np,nmem,que,qchem_Lsave):

    if Lusage:
        work = 'USAGE'
        name_class = classobj_dict[work]
        for key in name_class.__dict__.keys():
            print(f" {work}   \t:: {name_class.__dict__[key]}")
        return 0

    print("List this directory = ")
    mdir = os.path.dirname(__file__)        # __file__: this file location
    exe, mod, dirs, d_link = dir_all(mdir)
    sort_exe = sorted(exe)
    sort_mod = sorted(mod)
    sort_dir = sorted(dirs)
    #print(f"{exe} {mod}")
    if sort_dir:
        print("Directories:: ")
        for f in sort_dir:
            print(f"    {f}")

    if sort_exe:
        print("Executable:: ")
        if not work:
            for f in sort_exe:
                print(f"    {f}")
        else:
            for key in classobj_dict.keys():
                fkey = dir_classify(sort_exe, key, classobj_dict)
                if work:
                    name_class = classobj_dict[work]
                    for key in name_class.__dict__.keys():
                        if key in fkey:
                            print(f"    {key}.py\t:: {name_class.__dict__[key]}")
            print("  == remainder ")
            for f in sort_exe:
                print(f"    {f}")
    if sort_mod:       
        print("Module:: ")
        if not work:
            for f in sort_mod:
                print(f"    {f}")
        else:
            for key in classobj_dict.keys():
                fkey = dir_classify(sort_mod, key, classobj_dict)
                if work:
                    name_class = classobj_dict[work]
                    for key in name_class.__dict__.keys():
                        if key in fkey:
                            print(f"    {key}.py\t:: {name_class.__dict__[key]}")
            print("  == remainder ")
            for f in sort_mod:
                print(f"    {f}")
    if  work in lclass_var:
        if work == 'QChem':
            if fname == None:
                com = f"    qsub -N {quejob} -v qcjob=a.in -v nc={np} [-v opt=save] -pe numa {np} /gpfs/home/joonho/sandbox_gl/pypbs/sge_qchem.csh"
            else:
                com = f"    qsub -N {quejob} -v qcjob={fname} -v nc={np} "
                if qchem_Lsave:
                    com += '-v opt=save '
                com += f'-pe numa {np} /gpfs/home/joonho/sandbox_gl/pypbs/sge_qchem.csh'
            print(com)
            if not Lrun:
                print("\tuse -r to queue submit")
                print("    \"{}\" is aliased to \"{}\"".format(os.path.basename(__file__)+" -w QChem", "qchem_run"))
                print("    use :: qchem_run -j Queue-jobname -f qin-file -n num_process")
                print("    e.g.:: qchem_run -jopt -fCO2 -n4")
                print(" version:: use 3.1 for jmol and 5.1 for all the keywords")
                print("    For batch::\n\tdir_run.sh in \# $1 for qchem input type")
                if fname == None:
                    print(f"    in script::\n\t/gpfs/home/joonho/sciware/qchem5.1p/bin/qchem -np {np} a.in a.out")
                else:
                    print(f"    in script::\n\t/gpfs/home/joonho/sciware/qchem5.1p/bin/qchem -np {np} {fname}.in {fname}.out")
                print("\tenv variable such as $QC is not working in script")
            else:
                q = "will you run?"
                if yes_or_no(q):
                    os.system(com)
        elif work == 'AMP':
            print("SGE validation::")
            print(f"    qsub -v fname={fname} -v py_job=val -v scan=scan -v nc={np} -q {que} -pe numa {np} -l h_vmem={nmem} $SGE_HOME/sge_amp.csh")
        elif work == 'Qsleep':
            print("Dummy job::")
            print(f"            qsub -q {que} -pe numa {np} -l mem={nmem} $SB/pypbs/sge_sleep.csh")
            print("    Sample: qsub -q qname(sandy@slet02) -pe numa 6(num of process) -l mem=24G $SB/pypbs/sge_sleep.csh")
    if work in classobj_dict.keys():                    
        name_class = classobj_dict[work]
        for key in name_class.__dict__.keys():
            print(f" {work}   \t:: {name_class.__dict__[key]}")
    print(f"#Comment: -c    for classification'\
            \n\t  -u for usage: equal to -j USAGE\
            \n\t  -w {classobj_dict.keys()} for detail\
            ")

    return 0

def main():
    parser = argparse.ArgumentParser(description="display Usage for $SB/pypbs  ")
    parser.add_argument('-w','--work',  choices=['AMP','QChem','SGE','Qsub','Qsleep','GRMX'], help="several explanation option ")
    parser.add_argument('-f', '--infile', help='energy data input file')
    parser.add_argument('-r', '--run', action='store_true', help='run command')
    queue = parser.add_argument_group(title = 'SGE Queue')
    queue.add_argument('-j', '--que_job', default="opt", help='SGE Queue jobname')
    queue.add_argument('-q', '--queue', default="sandy@opt03", help='if you want to assign queue')
    queue.add_argument('-n', '--nprocess', default=12, type=int,  help='number of parallel process')
    queue.add_argument('-l', '--memory', default=2.0, type=float,  help='size of memory per process')
    qchem = parser.add_argument_group(title='Q-Chem')
    qchem.add_argument('-s', '--save', action='store_true', help='save option in Q-Chem run')
    parser.add_argument('-c','--classify', action='store_true', help="classify ")
    #parser.add_argument('-js','--specify', choices=['qcmo','nbo','eda'], help="present class details ")
    parser.add_argument('-u','--usage', action='store_true', help="present main details")
    args = parser.parse_args()

    works(args.work,args.usage,args.infile,args.run,args.que_job,args.nprocess,args.memory,args.queue,args.save)
    return 0

if __name__ == "__main__":
    main()
