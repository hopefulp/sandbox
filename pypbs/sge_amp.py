#!/home/joonho/anaconda3/bin/python

import argparse
import os
import sys
import numpy as np
from common import yes_or_no, list2str, whereami

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
def qsub_command(qjobname, ncore, mem, infile, amp_job, descriptor,p_function,p_minmax,nparam, Lte_force, str_hl, e_list, f_list, ntotal,ntr, ntype, nlist):
    if mem.isdigit():
        mem += 'G'
    elimit = e_list[0]
    comm  = f"qsub -N {qjobname}"
    comm += f" -pe numa {ncore} -l mem={mem} -v fname={infile} -v pyjob={amp_job} -v hl=\"{str_hl}\" -v el={elimit}"
    if Lte_force:
        comm += " -v tef=force"
    if descriptor:
        comm += f" -v des={descriptor}"
        comm += f" -v pf={p_function}"
        comm += f" -v pmm='{p_minmax[0]} {p_minmax[1]}'"
        comm += f" -v pn={nparam}"
    comm += " -v fl='" + ' '.join(f_list) + "'"
    comm += " -v nt=" + ' '.join(ntotal)
    comm += f" -v ntr={ntr}"
    comm += f" -v dtype={ntype}"
    comm += " -v dlist='" + ' '.join(nlist) + "'"
    comm += " $SB/pypbs/sge_amp.csh"
    #print(f"in qsub_command, {comm}")
    return comm

def submit(qjobname,ncore,infile,amp_job,descriptor,p_function,p_minmax,nparam,Lte_force,HL,elist,f_list,mem,ntotal,ntr,ntype,nlist):
    pwd = os.getcwd()
    str_hl = " ".join(str(x) for x in HL)
    #if yes_or_no("run ?"):
    #jobdir = make_dir(HL, elimit)
    elimit = elist[0]
    if elimit == 0.0:
        elimit=int(elimit)
    if float(f_list[0]) == 0.0:
        flimit=int(f_list[0])
    else:
        flimit=float(f_list[0])
    ### Make directory name!!
    jobdir = qjobname + 'E'+str(elimit)+'F'+str(flimit) # +'HL'+list2str(HL)
    #jobdir = 'NC' + str(ncore) + 'HL'+list2str(HL)+'E'+str(elimit)+'F'+str(flimit) 

    if jobdir in os.listdir(pwd):
        print(f"{jobdir} exists")
        new_dir = pwd + "/" + jobdir
        if yes_or_no("continue ?"):
            os.chdir(f"{new_dir}")
            ### goto comm
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
            ### goto comm
        else:
            com = "will you continue in here ?"
            if yes_or_no(com):
                print("not ready")
                sys.exit(11)
                #pass
            else:
                print("exit")
                sys.exit(10)
    """ Make dir for Test?
    if amp_job == 'te':
        os.mkdir('te')
        tedir = new_dir + "/" + "te"
        os.chdir(tedir)
        for amp_dir in amp_db:
            os.system(f"ln -s {pwd}/{amp_dir} {tedir}/{amp_dir}")
    """
    ###                   1        2     3      4       5       6           7           8       9       10    11   12     13    14
    #comm = qsub_command(qjobname, ncore, mem, infile, amp_job, descriptor, Lte_force, str_hl, elist, f_list, ntotal, ntr, ntype, nlist)
    comm = qsub_command(qjobname, ncore, mem, infile, amp_job, descriptor, p_function,p_minmax,nparam, Lte_force, str_hl, elist, f_list, ntotal, ntr, ntype, nlist)
    str1 = f"will you run? \n{comm}"
    if yes_or_no(str1):
        os.system(comm)
    ### Recover original directory
    os.chdir(f"{pwd}")
    return 0

def HL_list(hl_ini, hll):
    '''
    hl_ini: input HL
    hll[2]:    [total number of HL sets, step of number of nodes in each layer]
    '''
    hls = []
    for i in range(hll[0]):
        np_hl = np.array(hl_ini)
        np_hl += i * hll[1]      # *2
        hls.append(list(np_hl))
    print(f"{hls} in {whereami()}")
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

def amp_work(work,qjobname,ncore,infile,amp_job,descriptor,p_function,p_minmax,nparam,Lte_force,HL,hll,elist,f_list,mem,ntotal,ntr,ntype,nlist):
    ### Make subdir and run
    if work == 'scan':
        ### scan: multi HLs
        hl_lists = HL_list(HL, hll)
        for hlist in hl_lists:
            qname = get_qname_suff(qjobname, *hlist)
            submit(qname,ncore,infile,amp_job,descriptor,p_function,p_minmax,nparam,Lte_force,hlist,elist,f_list,mem,ntotal,ntr,ntype,nlist)
    elif work == 'onejob':
        qname = get_qname(qjobname, *HL)
        submit(qname,ncore,infile,amp_job,descriptor,p_function,p_minmax,nparam,Lte_force,HL,elist,f_list,mem,ntotal,ntr,ntype,nlist)
    ### make DATA BASE without making subdirectory
    else:   # work == 'db'
        str_hl = " ".join(str(x) for x in HL)
        ###                 1          2    3      4       5          6         7       8       9       10        11      12      13     14    15    16     17
        comm = qsub_command(qjobname, ncore, mem, infile, amp_job, descriptor,p_function,p_minmax,nparam, Lte_force, str_hl, elist, f_list, ntotal,ntr, ntype, nlist)
        ###
        str1 = f"will you run: \n{comm}"
        if yes_or_no(str1):
            os.system(comm)
    return 0

def main():
    parser = argparse.ArgumentParser(description="making new_dir & submit jobs in queue in SGE(mlet) ")
    parser.add_argument('-qj', '--queue_jobname', default='amp', help='job option:"train","test","md","validation","profile"')
    parser.add_argument('-m', '--mem', default='10G', help='memory usage')
    parser.add_argument('-i', '--input_file', default='OUTCAR', help='ASE readible file: extxyz, OUTCAR(VASP) ')
    parser.add_argument('-j', '--amp_job', default='tr', help='job option:"train","test","md","validation","profile"')
    ### Descriptor group
    descriptor_group = parser.add_argument_group(title="Descriptor generator")
    descriptor_group.add_argument('-des', '--descriptor', choices=['gs','zn','bs'], help="test new descriptor")
    descriptor_group.add_argument('-pf', '--param_function', default='log10', choices=['log10','powNN'], help="function for parameter interval")
    descriptor_group.add_argument('-pmm', '--param_minmax', nargs=2, default=[0.05, 5.0], help="min, max for param interval")
    descriptor_group.add_argument('-pn', '--nparam', default=4, help="num of parameters for descriptor")

    parser.add_argument('-tef', '--test_force', action='store_true', help='As for test: calculate force')
    parser.add_argument('-hl', '--hidden_layers', nargs='*', type=int, default=[8,8,8], help='Hidden Layer of lists of integer')
    w_group = parser.add_mutually_exclusive_group()
    w_group.add_argument('-s', '--scan', action='store_true', help='Scan Hidden Layer')
    w_group.add_argument('-db', '--ampdb', action='store_true', help='make amp_fingreprint.db')
    scan_group = parser.add_argument_group(title = 'control scan')
    scan_group.add_argument('-nhl', '--hl_number', type=int, default=5, help='total number of HL')
    scan_group.add_argument('-ihl', '--hl_inter', type=int, default=1, help='interval of HL')
    parser.add_argument('-nc', '--ncore', default=12, type=int, help='number of core needs to be defined')
    parser.add_argument('-el', '--e_convergence', nargs='+', default=['0.001','0.003'], help='energy convergence limit') # make default string
    parser.add_argument('-fl', '--f_convergence', nargs='*', default=['0.01', '0.04'],  help='force convergence limit, force coeff')
    data_group = parser.add_argument_group(title='Data')
    data_group.add_argument('-nt', '--ndata_total', nargs='*', help='cut total data: it requires two value for region')
    data_group.add_argument('-ntr', '--ndata_train', help='in case of te, define ND_train, ND_test is calculated')
    data_group.add_argument('-dtype','--data_s_type',default='pick',choices=['npart','int','div','pick'], help='data selection type: div-divide by dl[0] and remainder dl[1] for train, dl[2] for test ')
    data_group.add_argument('-dl','--data_s_list', nargs='+', help='data selection list')
    args = parser.parse_args()
    ### Amp work is ['scan','db','onejob']
    hll = []
    if args.scan:
        work = 'scan'
        hll = [args.hl_number, args.hl_inter]
    elif args.ampdb:
        work = 'db'
    else:
        work = 'onejob'
    ### run amp_work
            # 1             2              3             4             5              6               7                8                   9     
    amp_work(work,args.queue_jobname,args.ncore,args.input_file,args.amp_job,args.descriptor,args.param_function,args.param_minmax,args.nparam,args.test_force,args.hidden_layers,hll,args.e_convergence, args.f_convergence, args.mem, args.ndata_total, args.ndata_train, args.data_s_type, args.data_s_list)
    # 10                 11           12            13                  14             15             16               17               18     
    # 19
if __name__ == "__main__":
    main()
