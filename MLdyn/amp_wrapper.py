#!/home/joonho/anaconda3/bin/python

import numpy as np
import os
import re
import argparse
import sys
import time
import amp_ini
import amp_util
from ase.calculators.calculator import Parameters
from server_env import nXn
from common import whereami

def run_amp_wrapper(fin, job_submit, qname, ampjob, test_wrapper, key_values, hl, dlist, ndtotal, sftype, partition, nnode, nproc):
    '''
    sleeptime: There are two qsub(sbatch) for train/test. These have time interval using time.sleep()
    to test routine: control
        dlist = [1000, 1030, 3500, 3510]
        max_iter=10
        sleeptime = 10
    '''
    cwd = os.getcwd()
    cwdname = cwd.split('/')[-1]

    force_train = 1              ### None: Do not train force but test, 1/any to use default

    ### short running : -t r, running test can make database
    if test_wrapper == 'r':
        dlist       = [1000, 1030, 3500, 3510]
        max_iter    = 10
        sleeptime   = 10
        ncore       = 4
        mem         = '3G'
    elif test_wrapper == 'db':
        if not dlist:
            print("input data interval for making db")
            sys.exit(0)
        max_iter    = 10
        sleeptime   = 600
        if partition:
            ncore   =   nXn[partition]
        else:
            ncore   = 16
        mem         = '2.5G'
    else:
        ### main check point (1)
        if not dlist:
            dlist   = [1000,2000,3500,3800]      # [1000,1100,3500,4000] [1000,2500,3500,3600]
        max_iter    = 10000  # 5000(Nd300), 10000(Nd500, it takes one weak)
        sleeptime   = 600    # 600
        ### Nd500: pn 5 ~ 20: 
        ### Nd300: pn50 8.6G, in 96: np10 m8.6
        ncore       = '16'      # ncore 16, mem 11 for 188G
        mem         = '8G'     # 11G for ND500 15G

    if hl:
        str_hl = ' '.join(str(x) for x in hl)
    else:
        str_hl = '4'
        
    ### for more options: check amp_ini.Amp_string
    ### memory-ncore: for 188: 12G-12core for 96G: 12G-8core
    ### SET UP amp_run argument
    p_tr = Parameters({'jsubmit':job_submit, 'ampjob':'tr',  'ncore':ncore, 'hl':str_hl, 'elimit':'0.0001',
                                                            'mem':mem, 'dlist': dlist, 'queuename':qname+'tr',
                                                            'partition':partition, nnode:nnode})
    p_te = Parameters({'jsubmit':job_submit, 'ampjob':'te', 'ncore':'1', 'dlist': dlist, 'queuename':qname+'te',
                                                            'partition':partition, 'nnode':1})
    #print(f"p_tr partition {p_tr['partition']}")
    ### add keyword more
    if re.search('tr', ampjob):
        ### symmetry function needs to be modified in amp_ini.Amp_string.symmetry_function. if DB was made, this is not required
        ### sftype from args.sftyp so if not, it's None, max_iter is not input
        if sftype:
            p_tr['sftype'] = sftype
        else:
            p_tr['sftype'] = 'logm200n10'
        if 'max_iter' in locals():
            p_tr['max_iter'] = max_iter
        if nproc:
            p_tr['nproc'] = nproc
            print(f"nproc {p_tr.nproc} in {whereami()}")
    ### add keyword for tr and te
    if fin != 'OUTCAR':
        p_tr['fname'] = fin
        p_te['fname'] = fin
    if ndtotal:
        p_tr['ndtotal'] = ndtotal   # int
        p_te['ndtotal'] = ndtotal   # int
            
    if force_train == None:
        p_tr.train_f = None
        p_te.train_f = None
        
    ### make and run amp_run.py command for train
    if re.search('tr', ampjob):
        ### if amp_pot exist, do not run training
        if amp_util.is_amppot():    # does it need cwd?
            pass
        else:
            #if os.path.isfile("amp.amp") or os.path.isfile("amp-untrained-parameters.amp"):
            ### deprecate checking amp.pot
            #if amp_util.get_amppotname():
            #    print("There are amp potential file so exit")
            #    sys.exit(1)
            ### run training
            if job_submit == 'qsub':
                ampstr = amp_ini.Amp_string(add_amp_kw=p_tr)
                with open("mlet_tr.csh", 'w') as f:
                    f.write(ampstr.qscript)
            elif job_submit == 'sbatch':
                ampstr = amp_ini.Amp_string(add_amp_kw=p_tr)
                with open("sbatch_tr.sh", 'w') as f:
                    f.write(ampstr.qscript)
            elif job_submit == 'node':
                ampstr = amp_ini.Amp_string(add_amp_kw=p_tr)
            print(ampstr())
            ### check for show
            if test_wrapper == None or test_wrapper == 'db':
                os.system(ampstr())
            ### if check: print test and exit
            else:
                if re.search('te', ampjob):
                    ampstr = amp_ini.Amp_string(add_amp_kw=p_te)
                    if job_submit == 'qsub':
                        with open("mlet_te.csh", 'w') as f:
                            f.write(ampstr.qscript)
                    elif job_submit == 'sbatch':
                        with open("sbatch_te.sh", 'w') as f:
                            f.write(ampstr.qscript)
                    print(ampstr())
                sys.exit(2)
            ### wait after running amp_run.py
            while True:
                time.sleep(sleeptime)
                if not os.path.isfile("amp-log.txt"):
                    if job_submit == 'qsub':
                        print(f"{cwdname}-{qname}:: job is not loaded in queue")
                else:
                    print(f"{cwdname}-{qname}:: waiting until training finishes")
                ### this is not working: when calc.train() in amp_run.py fails, doesnot make amp_ini.amptrain_finish but lives untrained.amp
                if os.path.isfile("amp.amp") or os.path.isfile("amp-untrained-parameters.amp"):
                    break
            print("\n\n\ntraining is done")
    ############################################################################################    
    ### Test Part         
    if re.search('te', ampjob):
        if not amp_util.get_amppotname():
            print("Error:: No amp potfile so exit")
            sys.exit(10)
        if job_submit == 'qsub':
            ampstr = amp_ini.Amp_string(add_amp_kw=p_te)
            with open("mlet_te.csh", 'w') as f:
                f.write(ampstr.qscript)
        elif job_submit == 'node':
            ampstr = amp_ini.Amp_string(add_amp_kw=p_te)
        print(ampstr())
        ### check or show
        if test_wrapper == 'c' or test_wrapper == 's':
            sys.exit(3)
        else:
            if os.path.isfile(amp_ini.ampout_score):
                os.system(f"rm {amp_ini.ampout_score}") 
            os.system(ampstr())

        while True:
            time.sleep(sleeptime)
            print(f"{cwdname}-{qname}:: test is going until finding '{amp_ini.ampout_score}'")
            if os.path.isfile(amp_ini.ampout_score):
                break
        print("Job is Done")

    return 0

def main():
    parser = argparse.ArgumentParser(description='Process for Genetic Algorithm & Artificial Neural Network')
    parser.add_argument('-f', '--inf', default='OUTCAR', help='input file')
    parser.add_argument('-js', '--job_submit', default='qsub', choices=['qsub','node','sh','sbatch'], help='Select between node and login server')
    parser.add_argument('-qn', '--qname', default='test', help='jobname in queue')
    parser.add_argument('-j', '--job', default='trte', choices=['tr','te','trte'], help='ampjob')
    #parser.add_argument('-i', '--ini_setup', default = 'log10del12', choices=['log10del12'], help='input for amp_run.py')
    amp_gr = parser.add_argument_group(title='AMP')
    amp_gr.add_argument('-hl', '--hidden_layer', nargs='*', type=int, default=[4], help='Hidden Layer of lists of integer')
    amp_gr.add_argument('-dl', '--data_list', nargs='*', type=int, help='data index')
    amp_gr.add_argument('-nt', '--ndtotal', type=int, help='data index')

    amp_gr.add_argument('-sf', '--sym_function', help='input symmetry function with type of log10 and pow(N,N) such as logpn10max20, NN10')
    parser.add_argument('-t', '--test_wrapper', choices=['c','s','r','db'],  help='test amp_wrapper.py by "c" == "s" for check qsub file and r for short run w. qsub')
    parser.add_argument('-k', '--keyv', help='extra key values for special test')
    sbatch = parser.add_argument_group(title='SBATCH')
    sbatch.add_argument('-p', '--partition', default=1, type=int, help='designate partition for sbatch in iron')
    sbatch.add_argument('-nn', '--nnode', default=1, type=int, help='number of nodes')
    sbatch.add_argument('-np', '--nproc', default=0, type=int, help='number of nodes')

    args = parser.parse_args()

    run_amp_wrapper(args.inf, args.job_submit, args.qname, args.job, args.test_wrapper, args.keyv, args.hidden_layer, args.data_list, args.ndtotal, args.sym_function, args.partition, args.nnode, args.nproc)
    return 0


if __name__ == '__main__':
    main()

