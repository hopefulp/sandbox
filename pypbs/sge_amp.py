#!/home/joonho/anaconda3/bin/python

import argparse
import os
import sys
import numpy as np
from common import yes_or_no, list2str

amp_db = ['amp-fingerprint-primes.ampdb', 'amp-neighborlists.ampdb', 'amp-fingerprints.ampdb', 'OUTCAR']

def make_dir(HL, elimit):
    if HL:
        dname = "".join(str(HL))        # list to string
        print(dname)
        dname += str(elimit)
        return "HL"+dname
    else:
        print("HL is empty")
        sys.exit(10)
def amp_command(qjobname, ncore, mem, infile, amp_job, str_hl, elimit, f_list, ntotal, ntype, nlist):
    comm  = f"qsub -N {qjobname}"
    comm += f" -pe numa {ncore} -l mem={mem} -v fname={infile} -v pyjob={amp_job} -v hl=\"{str_hl}\" -v el={elimit}"
    comm += " -v fl='" + ' '.join(f_list) + "'"
    comm += f" -v nc={ncore}"
    comm += " -v ndata=" + ' '.join(ntotal)
    comm += f" -v dtype={ntype} "
    comm += " -v dlist='" + ' '.join(nlist) + "'"
    comm += " $SB/pypbs/sge_amp.csh"
    return comm

def submit(qjobname,ncore,infile,amp_job,nHL,elimit,f_list,mem,ntotal,ntr,ntype,nlist):
    pwd = os.getcwd()
    str_hl = " ".join(str(x) for x in nHL)
    #if yes_or_no("run ?"):
    #jobdir = make_dir(nHL, elimit)
    if elimit == 0.0:
        elimit=int(elimit)
    if float(f_list[0]) == 0.0:
        flimit=int(f_list[0])
    else:
        flimit=float(f_list[0])
    ### Make directory name!!
    jobdir = qjobname + 'E'+str(elimit)+'F'+str(flimit) # +'HL'+list2str(nHL)
    #jobdir = 'NC' + str(ncore) + 'HL'+list2str(nHL)+'E'+str(elimit)+'F'+str(flimit) 

    if jobdir in os.listdir(pwd):
        print(f"{jobdir} exists")
        if yes_or_no("continue ?"):
            pass
        else:
            sys.exit(1)
    else:
        com = f"will you make {jobdir}?"
        if yes_or_no(com):
            os.mkdir(jobdir)
            print(f"{jobdir} was made")
            new_dir = pwd + "/" + jobdir
            ### cd directory
            os.chdir(f"{new_dir}")
            for amp_dir in amp_db:
                os.system(f"ln -s {pwd}/{amp_dir} {new_dir}/{amp_dir}")
            #os.system(f"cp {pwd}/OUTCAR {new_dir}")
            print("amp- directories and OUTCAR were copied")
        else:
            print("exit")
            sys.exit(10)
        comm = amp_command(qjobname, ncore, mem, infile, amp_job, str_hl, elimit, f_list, ntotal, ntype, nlist)
        str1 = f"will you run: \n{comm}"
        if yes_or_no(str1):
            os.system(comm)
        ### Recover original directory
        os.chdir(f"{pwd}")
        return 0

def HL_list(nmax, *hl_ini):
    hls = []
    #nhl = nmax - int(hl_ini[0])
    nhl = nmax
    for i in range(nhl):
        np_hl = np.array(hl_ini)
        np_hl += i      # *2
        hls.append(list(np_hl))
    print(hls)
    return hls

def get_qname(queue_default, *hlist):
    if queue_default == 'amp':
        qname = queue_default+'HL'+list2str(hlist)
    else:
        qname = queue_default
    return qname
def get_qname_suff(queue_default, *hlist):
    qname = queue_default+'HL'+list2str(hlist)
    return qname

def amp_work(work, qjobname,ncore,infile,amp_job,nHL,max_hl,elimit,f_list,mem,ntotal,ntr,ntype,nlist):
    ### Make subdir and run
    if work == 'scan':
        hl_lists = HL_list(max_hl, *nHL)
        for hlist in hl_lists:
            qname = get_qname_suff(qjobname, *hlist)
            submit(qname,ncore,infile,amp_job,hlist,elimit,f_list,mem,ntotal,ntr,ntype,nlist)
    elif work == 'onejob':
        qname = get_qname(qjobname, *nHL)
        submit(qname,ncore,infile,amp_job,nHL,elimit,f_list,mem,ntotal,ntr,ntype,nlist)
    ### make database without making subdirectory
    else:
        str_hl = " ".join(str(x) for x in nHL)
        comm = amp_command(qjobname, ncore, mem, infile, amp_job, str_hl, elimit, f_list, ntotal, ntype, nlist)
        str1 = f"will you run: \n{comm}"
        if yes_or_no(str1):
            os.system(comm)
    return 0

def main():
    parser = argparse.ArgumentParser(description="submit jobs in queue in SGE(mlet) ")
    parser.add_argument('-qj', '--queue_jobname', default='amp', help='job option:"train","test","md","validation","profile"')
    parser.add_argument('-m', '--mem', default='10G', help='memory usage')
    parser.add_argument('-i', '--input_file', default='OUTCAR', help='ASE readible file: extxyz, OUTCAR(VASP) ')
    parser.add_argument('-j', '--amp_job', default='tr', help='job option:"train","test","md","validation","profile"')
    parser.add_argument('-hl', '--hidden_layers', nargs='*', type=int, default=[8,8,8], help='Hidden Layer of lists of integer')
    w_group = parser.add_mutually_exclusive_group()
    w_group.add_argument('-s', '--scan', action='store_true', help='Scan Hidden Layer')
    w_group.add_argument('-db', '--ampdb', action='store_true', help='make amp_fingreprint.db')
    parser.add_argument('-mh', '--max_hl', type=int, default=10, help='maximum number of HL')
    parser.add_argument('-nc', '--ncore', default=12, type=int, help='number of core needs to be defined')
    parser.add_argument('-el', '--e_convergence', default=0.001, type=float, help='energy convergence limit')
    parser.add_argument('-fl', '--f_convergence', nargs='*', default=[0.01, 0.01], help='force convergence limit')
    data_group = parser.add_argument_group(title='Data')
    data_group.add_argument('-nt', '--ndata_total', nargs='*', help='cut total data: it requires two value for region')
    data_group.add_argument('-ntr', '--ndata_train', help='in case of te, define ND_train, ND_test is calculated')
    data_group.add_argument('-dt','--data_s_type',default='pick',choices=['npart','int','div','pick'], help='data selection type: div-divide by dl[0] and remainder dl[1] for train, dl[2] for test ')
    data_group.add_argument('-dl','--data_s_list', nargs='+', help='data selection list')
    args = parser.parse_args()
    ### Amp work is ['scan','db','onejob']
    if args.scan:
        work = 'scan'
    elif args.ampdb:
        work = 'db'
    else:
        work = 'onejob'
    amp_work(work,args.queue_jobname,args.ncore,args.input_file,args.amp_job,args.hidden_layers,args.max_hl,args.e_convergence, args.f_convergence, args.mem, args.ndata_total, args.ndata_train, args.data_s_type, args.data_s_list)

if __name__ == "__main__":
    main()
