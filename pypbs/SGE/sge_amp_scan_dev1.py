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

def submit(qjobname,ncore,infile,amp_job,nHL,elimit,flimit,mem,ntotal,ntr,ntype,nlist):
    pwd = os.getcwd()
    str_hl = " ".join(str(x) for x in nHL)
    #if yes_or_no("run ?"):
    #jobdir = make_dir(nHL, elimit)
    if elimit == 0.0:
        elimit=int(elimit)
    if flimit == 0.0:
        flimit=int(flimit)
    ### Make directory name!!
    #jobdir = 'HL'+list2str(nHL)+'E'+str(elimit)+'F'+str(flimit)
    jobdir = 'NC' + str(ncore) + 'HL'+list2str(nHL)+'E'+str(elimit)+'F'+str(flimit) 

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
        comm  = f"qsub -N {qjobname}"
        comm += f" -pe numa {ncore} -l mem={mem} -v fname={infile} -v pyjob={amp_job} -v hl=\"{str_hl}\" -v el={elimit} -v fl={flimit} -v nc={ncore}"
        comm += " -v ndata=" + ' '.join(ntotal)
        comm += f" -v dtype={ntype} "
        comm += " -v dlist='" + ' '.join(nlist) + "'"
        comm += " $SB/pypbs/sge_amp.csh"
        str1 = f"will you run: \n{comm}"
        if yes_or_no(str1):
            os.system(comm)
        ### Recover original directory
        os.chdir(f"{pwd}")
        return 0
def HL_list(nmax, *hl_ini):
    hls = []
    nhl = nmax - int(hl_ini[0])
    for i in range(nhl+1):
        np_hl = np.array(hl_ini)
        np_hl += i
        hls.append(list(np_hl))
    print(hls)
    return hls

def get_qname(queue_default, hlist):
    if queue_default == 'amp':
        qname = queue_default+'HL'+list2str(hlist)
    else:
        qname = queue_default
    return qname
def get_qname_app(queue_default, hlist):
    qname = queue_default+'HL'+list2str(hlist)
    return qname

def main():
    parser = argparse.ArgumentParser(description="submit jobs in queue in SGE(mlet) ")
    parser.add_argument('-qj', '--queue_jobname', default='amp', help='job option:"train","test","md","validation","profile"')
    parser.add_argument('-m', '--mem', default='10G', help='memory usage')
    parser.add_argument('-f', '--input_file', default='OUTCAR', help='ASE readible file: extxyz, OUTCAR(VASP) ')
    parser.add_argument('-j', '--amp_job', default='tr', help='job option:"train","test","md","validation","profile"')
    parser.add_argument('-hl', '--hidden_layers', nargs='*', type=int, default=[8,8,8], help='Hidden Layer of lists of integer')
    w_group = parser.add_mutually_exclusive_group()
    w_group.add_argument('-s', '--scan', action='store_true', help='Scan Hidden Layer')
    w_group.add_argument('-db', '--ampdb', action='store_true', help='make amp_fingreprint.db')
    parser.add_argument('-mh', '--max_hl', type=int, default=10, help='maximum number of HL')
    parser.add_argument('-nc', '--ncore', default=12, type=int, help='number of core needs to be defined')
    parser.add_argument('-el', '--e_convergence', default=0.001, type=float, help='energy convergence limit')
    parser.add_argument('-fl', '--f_convergence', default=0.01, type=float, help='force convergence limit')
    data_group = parser.add_argument_group(title='Data')
    data_group.add_argument('-nt', '--ndata_total', nargs='*', help='cut total data: it requires two value for region')
    data_group.add_argument('-ntr', '--ndata_train', help='in case of te, define ND_train, ND_test is calculated')
    data_group.add_argument('-dt','--data_s_type',default='pick',choices=['npart','int','div','pick'], help='data selection type: div-divide by dl[0] and remainder dl[1] for train, dl[2] for test ')
    data_group.add_argument('-dl','--data_s_list', nargs='+', help='data selection list')
    args = parser.parse_args()
    
    if args.scan:
        work = 'scan'
    elif args.ampdb:
        work = 'db':
    else:
        work = 'onejob'
    hl_ini=args.hidden_layers
    if args.scan:
        hl_lists = HL_list(args.max_hl, *hl_ini)
        for hlist in hl_lists:
            qname = get_qname_list(args.queue_jobname, hlist)
            submit(qname, args.ncore, args.input_file, args.amp_job, hlist, args.e_convergence, args.f_convergence, args.mem, args.ndata_total, args.ndata_train, args.data_s_type, args.data_s_list)
    else:
        qname = get_qname(args.queue_jobname, hl_ini)
        submit(qname, args.ncore, args.input_file, args.amp_job, hl_ini, args.e_convergence, args.f_convergence, args.mem, args.ndata_total, args.ndata_train, args.data_s_type, args.data_s_list)

if __name__ == "__main__":
    main()
