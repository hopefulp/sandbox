#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os
import socket
from common import *

def print_sge(inf,software,qjobname,sub_job,np,mem,Lscan_saver,data_int,hl,el,fl):
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
        print("Usage:: qsub -N qjob -v qcjob=infile -v ver='3.2p' -l mem=3G /gpfs/home/joonho/sandbox_gl/pypbs/sge_qchem.csh")
        print("        -N qjob: queue job name -j")
        print("        -v qcjob: qchem input file, use -i")
        print("        -v ver: version ['3.2p','4.3s', '5.1p'] parallel and serial")
        print("        -v save=ok: for save option for qchem")
        q = "Do you want to see sge_qchem script?"
        if not Lscan_saver:
            com = f"qsub -N {qjobname} -pe numa {np} -l mem={mem}G -v qcjob={inf} -v np={np} $SB/pypbs/sge_qchem.csh"
        else:
            com = f"qsub -N {qjobname} -pe numa {np} -l mem=2G -v qcjob={inf} -v np={np} -v scan=ok $SB/pypbs/sge_qchem.csh"
    elif software == 'vasp':
        if not qjobname:
            print("qsub -N pe500 -pe numa %d -v np=%d -v dir=pe500 $SB/pypbs/sge_vasp.csh" % (np))
            print("use -d qjobname -n np")
        else:
            com = "qsub -N %s -pe numa %d -v np=%d -v dir=%s $SB/pypbs/sge_vasp.csh" % (qjobname,np,np,qjobname)
            print(com)
    elif software == 'sleep':
        if qjobname:
            com = "qsub -N %s -pe numa %d $SB/pypbs/sge_sleep.csh" % (qjobname, np)
        else:
            com = "qsub -pe numa %d $SB/pypbs/sge_sleep.csh" % (np)
        print(com)

    elif software == 'amp':
        com = f"qsub -N {qjobname} -pe numa {np} -v fname=OUTCAR -v np={np} -v pyjob={sub_job} "
        if mem > 2.0:
            com += f"-l mem={mem}G "
        if data_int:
            st = " ".join(data_int)
            com += f"-v di=\"{st}\" "
        if hl:
            st = " ".join(hl)
            com += f"-v hl=\"{st}\" "
        if el:
            com += f"-v el={el} "
        if fl:
            com += f"-v fl={fl} "
        if Lscan_saver:
            com += "-v scan=ok "
        com += "$SB/pypbs/sge_amp.csh"
        #print(com) 
    if 'com' in locals():
        print(com)
        if yes_or_no("would you want to run?"):
            os.system(com)
    return 0

def print_chi(inf, software, np, Lscan_saver):
    _HOME = os.getenv('HOME')

    if software == 'qchem':
        print("Check for library: ")
        com = "module li"
        os.system(com)
        print("Serial::\n\t$ qchem a.in a.out")
        print("Parallel::\n\t$ mpirun -np n_process $QC/exe/qcprog.exe a.in $QCSCRATCH/{inf} > a.out &")

        com = f"mpirun -np {np} $QC/exe/qcprog.exe {inf}.in $QCSCRATCH/{inf} > {inf}.out &"
        if yes_or_no("Do you want to run:\n " + com):
            os.system(com)
    return 0

def job_description(fname, software, qjobname, sub_job, np, mem, scan_saver, data_int, hl, el, fl):
    server =  socket.gethostname()
    if server == 'login':
        print_sge(fname, software, qjobname, sub_job, np, mem, scan_saver,data_int,hl,el,fl) 
    elif server == 'chi':
        print_chi(fname, software, np, scan_saver)
    else:
        print(f"hostname error {server}")
        sys.exit(1)
    return 0        

def main():
    parser = argparse.ArgumentParser(description="how to use qsub:\n Usage:: qsub_server.py sge qchem -j qjobname(dirname) -i file.in -n np")
    #parser.add_argument('-s'. '--server', default='sge', choices=['sge', 'chi'], help='jobname in pbs file')
    parser.add_argument('software', choices=['qchem', 'grmx', 'vasp','sleep','amp'], help='kind of software')
    parser.add_argument('-qj', '--qjobname', help='qjob name for qsub|directory name')
    parser.add_argument('-js', '--sub_job', default='tr', help='for amp, [tr, te]')
    parser.add_argument('-i', '--inf', help='input filename')
    parser.add_argument('-m', '--mem', default=2.0, type=float, help='memory size for qsub')
    parser.add_argument('-n', '--np', default=16, type=int, help='number of process')
    gr_amp = parser.add_argument_group()
    gr_amp.add_argument('-sc', '--scan', action='store_true', help='scan amp for several Hidden Layer|save for qchem')
    gr_amp.add_argument('-di', '--data_interval', nargs="*", help='data selection for training and test')
    gr_amp.add_argument('-hl', '--hidden_layer', nargs="*", help='hidden layer')
    gr_amp.add_argument('-el', '--e_limit', type=float, help='training energy accuracy')
    gr_amp.add_argument('-fl', '--f_limit', type=float, help='training force accuracy')

    args = parser.parse_args()
    
    if args.inf.endswith('.in'):
        fname = args.inf[:-3]
    else:
        fname = args.inf

    job_description(fname, args.software, args.qjobname, args.sub_job, args.np, args.mem, args.scan, args.data_interval, args.hidden_layer, args.e_limit,args.f_limit) 

if __name__ == '__main__':
    main()
