#!/home/joonho/anaconda3/bin/python

import numpy as np
import os
import argparse
import sys
import time
import amp_ini
from ase.calculators.calculator import Parameters



def run_amp_wrapper(fin, job_submit, Lcheck_string, Lrun_test, hl, qname):
    '''
    There are two qsub for train/test. These have time interval using time.sleep()
        import time,
    not yet: sleep is not loaded in sge-script.csh using cshell
    '''
    cwd = os.getcwd().split('/')[-1]
    ### main check point (1)
    force_train = 1              ### None: Do not train force but test, 1/any to use default
    sleeptime = 300
    dlist = [1000,2500,3500,3600]      # [1000,1100,3500,4000]
    if hl:
        str_hl = ' '.join(str(x) for x in hl)
    else:
        str_hl = '4'
        
    ### for more options: check amp_ini.Amp_string

    ### SET UP amp_run argument
    p1 = Parameters({'jsubmit':job_submit, 'ampjob':'tr',  'ncore':'14', 'max_iter':2000, 'hl':str_hl, 'elimit':'0.0001',
                    'dlist':dlist, 'mem':'12G', 'queuename':qname+'tr'})
    p2 = Parameters({'jsubmit':job_submit, 'ampjob':'te', 'ncore':'1', 'dlist':dlist, 'queuename':qname+'te'})

    if force_train == None:
        p1.train_f = None
        p2.train_f = None
        
    ### make and run amp_run.py command for train
    if os.path.isfile("amp.amp") or os.path.isfile("amp-untrained-parameters.amp"):
        print("There are amp potential file so exit")
        sys.exit(1)
    else:
        if job_submit == 'qsub':
            ampstr = amp_ini.Amp_string(add_amp_kw=p1)
            with open("mlet_tr.csh", 'w') as f:
                f.write(ampstr.qscript)
            str_2run = ampstr()
            
        elif job_submit == 'node':
            ampstr = amp_ini.Amp_string(add_amp_kw=p1)
            str_2run = "amp_run.py" + ampstr() + " -g &"
        print(str_2run)
        if not Lrun_test and not Lcheck_string:
            os.system(str_2run)
        ### if check: print test and exit
        elif Lcheck_string:
            ampstr = amp_ini.Amp_string(add_amp_kw=p2)
            if job_submit == 'qsub':
                with open("mlet_te.csh", 'w') as f:
                    f.write(ampstr.qscript)
                str_2run = ampstr()
            else:
                str_2run = "amp_run.py" + ampstr() + " -g &"
            print(str_2run)
        if Lcheck_string: sys.exit(2)
        ### make and run amp_run.py command for test
        while True:
            time.sleep(sleeptime)
            if not os.path.isfile("amp-log.txt"):
                if job_submit == 'qsub':
                    print(f"{cwd}:: job is not loaded in queue")
            else:
                print(f"{cwd}:: waiting until training finishes")
            if os.path.isfile("amp.amp") or os.path.isfile("amp-untrained-parameters.amp"):
                break
        print("\n\n\ntraining is done")
    ############################################################################################    
    ### Test Part                
    if job_submit == 'qsub':
        ampstr = amp_ini.Amp_string(add_amp_kw=p2)
        with open("mlet_te.csh", 'w') as f:
            f.write(ampstr.qscript)
        str_2run = ampstr()

    elif job_submit == 'node':
        ampstr = amp_ini.Amp_string(add_amp_kw=p2)
        str_2run = "amp_run.py" + ampstr() + " -g &"
    print(str_2run)
    if os.path.isfile(amp_ini.ampout_te_f_chk):
        os.system(f"rm {amp_ini.ampout_te_f_chk}") 
    os.system(str_2run)

    while True:
        time.sleep(sleeptime)
        print(f"{cwd}:: test is going until finding '{amp_ini.ampout_te_f_chk}'")
        if os.path.isfile(amp_ini.ampout_te_f_chk):
            break
    print("Job is Done")

    return 0

def main():
    parser = argparse.ArgumentParser(description='Process for Genetic Algorithm & Artificial Neural Network')
    parser.add_argument('-f', '--inf', default='OUTCAR', help='input file')
    parser.add_argument('-js', '--job_submit', default='qsub', choices=['qsub','node','sh'], help='Select between node and login server')
    #parser.add_argument('-i', '--ini_setup', default = 'log10del12', choices=['log10del12'], help='input for amp_run.py')
    parser.add_argument('-hl', '--hidden_layer', nargs='*', type=int, default=[4], help='Hidden Layer of lists of integer')
    parser.add_argument('-qn', '--qname', default='test', help='jobname in queue')
    parser.add_argument('-c', '--check_string', action='store_true', help='check command of amp_run.py but not submit the jobs')
    parser.add_argument('-te', '--test_only', action='store_true', help='check command of amp_run.py but not submit the jobs')

    args = parser.parse_args()

    run_amp_wrapper(args.inf, args.job_submit, args.check_string, args.test_only, args.hidden_layer, args.qname)
    return 0


if __name__ == '__main__':
    main()

