#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os
from common import *

def print_sge(software,dname,jobname,np,Lscan_saver):
    _HOME = os.getenv('HOME')
    if not software:
        print("Use:: qfree - to check freed node")
        print("      qhist - to check user occupancy")
        print("      qmem  - to check memory usage")
        print("      qstat - see my job")
        print("            -f for full")
        print("            -f | sed '/---/d' | sed 's*0/**' | sed 's*resv/**' ")
        print("Usage:: -s software to see how to run qsub {qchem, grmx, vasp, sleep, amp}")

    elif software == 'qchem':
        print("SGE(Rhee's group): Q-Chem")
        print("Usage:: qsub -N jobname -v qcjob=qchem_input -v ver='3.2p' /gpfs/home/joonho/sandbox_gl/pypbs/sge_qchem.tcsh")
        print("        -N jobname: listed in 'qstat', jobname should be alphanumeric, use -d")
        print("        -v job: variable in 'qchem job.in job.out' in script sge_qchem.csh, use -j")
        print("        -v ver: version ['3.2p','4.3s', '5.1p'] parallel and serial")
        print("        -v save=ok: for save option for qchem")
        q = "Do you want to see sge_qchem script?"
        #if yes_or_no(q):
        #    sandbox = _HOME + '/sandbox_gl/pypbs/sge_qchem.csh'
        #    com = "more %s" % sandbox
        #    os.system(com)
        if not Lscan_saver:
            com = f"qsub -N {dname} -pe numa {np} -v qcjob={jobname} -v np={np} $SB/pypbs/sge_qchem.csh"
        else:
            com = f"qsub -N {dname} -pe numa {np} -v qcjob={jobname} -v np={np} -v scan=ok $SB/pypbs/sge_qchem.csh"

    elif software == 'vasp':
        if not dname:
            print("qsub -N pe500 -pe numa %d -v np=%d -v dir=pe500 $SB/pypbs/sge_vasp.csh" % (np))
            print("use -d dname -n np")
        else:
            com = "qsub -N %s -pe numa %d -v np=%d -v dir=%s $SB/pypbs/sge_vasp.csh" % (dname,np,np,dname)
            print(com)
    elif software == 'sleep':
        if dname:
            com = "qsub -N %s -pe numa %d $SB/pypbs/sge_sleep.csh" % (dname, np)
        else:
            com = "qsub -pe numa %d $SB/pypbs/sge_sleep.csh" % (np)
        print(com)
    elif software == 'amp':
        if Lscan_saver:
            com = f"qsub -N {dname} -pe numa {np} -v fname=OUTCAR -v np={np} -v pyjob=tr $SB/pypbs/sge_amp.csh"
        else:
            com = f"qsub -N {dname} -pe numa {np} -v fname=OUTCAR -v np={np} -v pyjob=tr $SB/pypbs/sge_amp.csh"
        print(com) 
    if 'com' in locals():
        print(com)
        if yes_or_no("would you want to run?"):
            os.system(com)
    return 0

def print_chi(software):
    _HOME = os.getenv('HOME')
    if software == 'qchem':
        print("Check for library: ")
        com = "module li"
        os.system(com)
        print("If serial   is loaded: $qchem a.in a.out")
        print("If parallel is loaded: $mpirun -np n_process $QC/exe/qcprog a.in $QCSCRATCH > a.out &")
    return 0

def job_description(server, software, dname, jname, np, amp_Lscan):
    if server == 'sge':
        print_sge(software,dname,jname, np, amp_Lscan) 
    elif server == 'chi':
        print_chi(software)
    return 0        

def main():
    parser = argparse.ArgumentParser(description='how to use qsub')
    parser.add_argument('server', default='sge', nargs='?', choices=['sge', 'chi'], help='jobname in pbs file')
    parser.add_argument('software', choices=['qchem', 'grmx', 'vasp','sleep','amp'], help='kind of software')
    parser.add_argument('-d', '--dirname', help='job directory name|job name for qsub')
    parser.add_argument('-j', '--jobname', help='additional job name for qchem')
    parser.add_argument('-n', '--np', default=16, type=int, help='number of process')
    gr_amp = parser.add_argument_group()
    gr_amp.add_argument('-sc', '--scan', action='store_true', help='scan amp for several Hidden Layer|save for qchem')
    args = parser.parse_args()

    job_description(args.server, args.software, args.dirname, args.jobname, args.np, args.scan) 

if __name__ == '__main__':
    main()
