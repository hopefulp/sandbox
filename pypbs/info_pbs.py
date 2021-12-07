#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_all, dir_classify_n, whereami, MyClass
#from common import MyClass_t as MyClass
#from varname import nameof

### thest are globals() not locals()|vars()
sge     = MyClass('sge')
qsub    = MyClass('qsub')
ssh     = MyClass('ssh')
grmx    = MyClass('grmx')
usage   = MyClass('usage')
amp     = MyClass('amp')
qchem   = MyClass('qchem')
qsleep  = MyClass('qsleep')
process = MyClass('process')
slurm   = MyClass('slurm')
kisti   = MyClass('kisti')

# if not use input variables, define here
sge.usage="$qstat -f    # see all nodes(/used/total) and my job\
            \n\t\t$qstatf      # aliast qstat -f ...\
            \n\t\t$qfree       # free nodes\
            \n\t\t$qfree       # free nodes\
            \n\t\t$qhist       # shows number of running cores and jobs in Queue\
            \n\t\t$qmem        # memory\
            "
sge.sge_dir_run     = "\n\tsge_dir_run.py\
                    \n\t    run in a directory w/wo qsub\
                    "
qsub.qsub_server     = f"QSUB JOB\
                    \n\tautodetection of server ['login','chi']\
                    \n\tWorks: ['qchem','amp','vasp','sleep']\
                    \n\t    Q-Chem:\
                    \n\t\tOne Job::\
                    \n\t\t    qsub_server.py qchem -qj Nimono -i CO2M06(.in) -n 16 -m 5(G)\
                    \n\t\t    qsub_server.py qchem -qj CC6 -i 6-CC-NiFe.in -n 6 -m 6 -no skylake@node14\
                    \n\t\tN.B.::\
                    \n\t\t    showall.py -s -j qchem\
                    \n\t\t    job downs in opt&slet -> specify -no skylake@node??\
                    \n\t\t\tby checking qmem & qstat\
                    \n\tOptions:\
                    \n\t    PBS:\
                    \n\t\t-qj --qjobname\
                    \n\t\t-i --inf inputfile\
                    \n\t\t-m --mem memory\
                    \n\t\t-n --np  nproc\
                    \n\t\t-no --node specify node\
                    "
qsub.usage          ="qsub -N jobname -v var=a_variable /gpfs/home/joonho/sandbox_gl/pypbs/sge_qchem.csh\
                    \n\t\t$var can be used as variable in the script\
                    \n\t\tscript name should be full name\
                    "
qsub.qdel           ="qdel (from) job-ID (to) job-ID\
                    \n\tto kill process in master node: rt. process.prockill\
                    "
process.prockill    ="kill (from) PID (to) PID\
                    "
ssh.ssh_node_byproc = "\n\tssh_node_byproc.sh node index_proc num_proc[default=16]\
                    \n\t    kill process in a node\
                    "
ssh.ssh_mlet_scan_nodes   = "\n\tssh_mlet_nodes.sh [ls|rm|mkdir|qchem|ln|ps|psl] [$2]\
                    \n\t    : scan all nodes in mlet\
                    \n\t    depending on $1, $2 is required\
                    \n\tUsage::\
                    \n\t    psl python - runs ssh all-nodes ps aux | grep python\
                    \n\t\t:ps shows all the line, psl shows just the number of lines for python process\
                    \n\t:: also check alias.sh for chechnodes for qstat nodes and checknode for one node\
                    "
ssh.ssh_node_kill   = "\n\tssh_node_kill.sh proc_name nodes\
                    \n\t    : for the given nodes find process of proc_name the kill each process\
                    "
grmx.usage="qsub -v tpr=mdname sge_mdrun2.sh"
usage.gen="#Comment for sge: -w [sge|amp|grmx|qchem|sleep|qsub for arguments]\
            \n\t\tAMP::     sge_amp.tcsh\
            \n\t\tGromacs:: sge_mdrun2.tcsh\
            \n\t\tQ-Chem::  sge_qchem.csh\
            \n\t\tdummy job::  sge_sleep.csh"
usage.files = ['u.py']
amp.usage="\n\t\tpy_job for [validation||train||test]\
            \n\t\tscan variable (if exists) average for several trial with the same [HL, ene convergency]\
            \n\t\tsetenv SGE_HOME $HOME/sandbox_gl/pypbs\
            \n\t\tfind node:: sge_nodes.py || qstat -f\
            \n\tin script::\n\tsetenv PYTHONPATH $HOME/sandbox_gl/pycommon:$HOME/sandbox_gl/mymplot\
            \n\t\tset PYTHON = \"$HOME/anaconda3/bin/python\"\
            \n\t\tset EXE = \"$HOME/sandbox_gl/py_ai/amp_ene.py\"\
            \n\t\t$PYTHON $EXE $fname $py_job -hl 4 4 4 -el 0.001 -n 5 -g\
            "
amp.sge_amp="sge_amp.py \
             sge_amp.csh\
             "
amp.amp_scan="amp_scan.sh\
            \n\t\tdirect run at node with scan\
            "
# kisti  = " " makes error: 149: gkey==instance.name if not instance
kisti.pbs_vasp_kisti_skl = "Usage:: python $sbvas/vasp_script.py args\
                            \n\t\tSkylake for VASP\
                            \n\t\tsetting: nnode=20, nproc=40 in 40\
                            " 
kisti.pbs_vasp_kisti_skl2 = "setting: nnode=20, nproc=20 in 40 for memory overflow\
                            "  

kisti.pbs_vasp_kisti_knl = "setting: nnode=10, nproc=68 in 68\
                            \n\t\tKnight for parallel process such as LAMMPS\
                            "  
kisti.pbs_vasp_kisti_knl2 = "setting: nnode=10, nproc=64 in 68\
                            \n\t\tfor fast calculation?\
                            "  


qchem.usage=""
qsleep.usage=""

slurm.slurm_sbatch      ="link file to slurm_sbatch_vasp.sh"
slurm.slurm_sbatch_vasp ="to run vasp in slurm\
                        \n\tUse: command.py slurm -sj vasp\
                        "
slurm.slurm_sbatch_py   ="to run python in node\
                        \n\tUse: command.py slurm -sj mldyn\
                        "

classobj_dict={'SGE': sge, 'GRMX': grmx, 'AMP':amp, 'USAGE': usage, 'Qsub': qsub, 'QChem': qchem, 'Qsleep': qsleep}
lclass_var=['QChem', 'AMP', 'Qsleep']
#print(classobj_dict['SGE'].usage)

def works(Lclassify, work,Lusage,fname,Lrun,quejob,np,nmem,que,qchem_Lsave):

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
        if not Lclassify:
            for f in sort_exe:
                print(f"    {f}")
        ### classify
        else:
            ### instance.name is string, stored as self.name in MyClass
            ### class instance variable (sge) appears in the same letter as instance.name='sge'
            #print(MyClass.instances)
            for instance in MyClass.instances:
                #print(globals().keys())
                for gkey in globals().keys():
                    if gkey == instance.name:    # 'sge'(class instance) == 'sge'(string)
                        break
                if work != instance.name:
                    ckeys = dir_classify_n(sort_exe, instance.name, globals()[gkey], Lwrite=1)
                else:
                    ckeys = dir_classify_n(sort_exe, instance.name, globals()[gkey], Lwrite=0)
                    for ckey in ckeys:
                        print(f"    {ckey}.py[sh]\t:: {globals()[gkey].__dict__[ckey]}")
            print("  Remainder == ")
            for f in sort_exe:
                print(f"    {f}")
    if sort_mod:       
        print("Module:: ")
        if not Lclassify:
            for f in sort_mod:
                print(f"    {f}")
        else:
            for instance in MyClass.instances:
                for gkey in globals().keys():
                    if gkey == instance.name:
                        break
                if not work or work != instance.name:
                    ckeys = dir_classify_n(sort_mod, instance.name, globals()[gkey], Lwrite=1)
                else:
                    ckeys = dir_classify_n(sort_mod, instance.name, globals()[gkey], Lwrite=0)
                    for ckey in ckeys:
                        print(f"    {ckey}.py\t:: {globals()[gkey].__dict__[ckey]}")
            print("  Remainder == ")
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
    #if work in classobj_dict.keys():
    print("Instances:: ", end='')
    for instance in MyClass.instances:
        print(f"{instance.name}", end=' ') 
    print("\n\t    -w for detail")

    return 0

def main():
    parser = argparse.ArgumentParser(description="PBS in $SB/pypbs  ")
    parser.add_argument('-c','--classify', action='store_false', help="classify ")
    #parser.add_argument('-w','--work',  choices=['AMP','QChem','SGE','Qsub','Qsleep','GRMX'], help="several explanation option ")
    parser.add_argument('-w','--work',  help="several explanation option ")
    parser.add_argument('-f', '--infile', help='energy data input file')
    parser.add_argument('-r', '--run', action='store_true', help='run command')
    queue = parser.add_argument_group(title = 'SGE Queue')
    queue.add_argument('-j', '--que_job', default="opt", help='SGE Queue jobname')
    queue.add_argument('-q', '--queue', default="sandy@opt03", help='if you want to assign queue')
    queue.add_argument('-n', '--nprocess', default=12, type=int,  help='number of parallel process')
    queue.add_argument('-l', '--memory', default=2.0, type=float,  help='size of memory per process')
    qchem = parser.add_argument_group(title='Q-Chem')
    qchem.add_argument('-s', '--save', action='store_true', help='save option in Q-Chem run')
    #parser.add_argument('-js','--specify', choices=['qcmo','nbo','eda'], help="present class details ")
    parser.add_argument('-u','--usage', action='store_true', help="present main details")
    args = parser.parse_args()

    works(args.classify,args.work,args.usage,args.infile,args.run,args.que_job,args.nprocess,args.memory,args.queue,args.save)
    return 0

if __name__ == "__main__":
    main()
