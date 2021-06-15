#!/home/joonho/anaconda3/bin/python

import numpy as np
import os
from ga import GA
import glob
import argparse
#from gaamp import qsub_restart
import sys
from common import whereami 
from ase.calculators.calculator import Parameters

def calc_ga(job_submit, exe, fin, fout,ndata, max_iter, setup, nchromo, nhl, nnode, ngenerations, nparent_sets, mutation_percent, Lcheck_string, Ltest):

    calc = GA(setup=setup, nchromo=nchromo, nhl=nhl, nnode=nnode)
    ### whether train force or not
    force_train = 1         # [1,any|None], True: force_train ; None: energy-only train

    ### run_generation params
    start_point = 0
    subdir_prefix = 'ch'

    if Ltest or Lcheck_string:
        sleeptime = 10
        p = Parameters({'jsubmit':job_submit, 'max_iter':300 })
    else:
        sleeptime = 300
        p = Parameters({'jsubmit':job_submit,'max_iter':max_iter })
    ### include force or not
    if force_train == None:     # else default p.train_f is set
        p.train_f = None
    if  ndata:
        p.ndtotal = ndata
    print(f"ngenerations {ngenerations}")
    # make jobstring and pass to calc or pass parameters then calc.run_generation make the jobstring due to hl made in run_generation
    calc.run_generation(job_submit=job_submit, dir_prefix=subdir_prefix, istart=start_point,
                        ngenerations=ngenerations, nparent_sets=nparent_sets, mutation_percent=mutation_percent, 
                        check_print=Lcheck_string, sleeptime=sleeptime,
                        ml_job_setting=p)

    print('GA CALCULATION IS DONE')
    return 0

def main():
    parser = argparse.ArgumentParser(description='Process for Genetic Algorithm & Artificial Neural Network')
    parser.add_argument('-js', '--job_submit', default='node', choices=['qsub','node'], help='Select between node and login server')
    parser.add_argument('-e', '--exe', default='mllorenz.py', help='python script for ML')
    parser.add_argument('-if', '--ifile', default='OUTCAR', help='input file')
    parser.add_argument('-of', '--ofile', help='save pt model file')
    parser.add_argument('-u', '--usage', action='store_true', help='show usage')
    parser.add_argument('-c', '--check_string', action='store_true', help='check command string of amp_run.py but not submit the jobs')
    parser.add_argument('-t', '--Ltest', action='store_true', help='run test job for checking source code')
    genetic_params = parser.add_argument_group(title = 'Genetic Algorithm')
    genetic_params.add_argument('-i', '--ini_setup', default = 'random', choices=['random','reference','ongoing'], help='reference/ongoing from GANNtest file')
    genetic_params.add_argument('-nch', '--nchromo', default=20, type=int, help="number of HLs in each generation")
    genetic_params.add_argument('-hl', '--hlayers', default=10, type=int, help="fixed number of layers in each HL")
    genetic_params.add_argument('-nn', '--nnodes', default=20, type=int, help="number of nodes in each HL")
    genetic_params.add_argument('-mi', '--max_iter', default=1000000, type=int, help="number of maximum iteration")
    genetic_params.add_argument('-nd', '--ndata', type=int, default=80000,  help='total number of data')
    genetic_params.add_argument('-np', '--nparents_sets', default=4, type=int, help="number of parents mating: select best hls")
    genetic_params.add_argument('-mp', '--prob_mutation', default=10, type=int, help="percentage of mutation")
    genetic_params.add_argument('-ng', '--ngenerations', default=10, type=int, help="number of generation")

    args = parser.parse_args()
    if args.usage:
        print("For Node")
        print(f"{os.path.basename(__file__)} -e 'mllorenz.py' -js node -nch 20 -hl 10 -nn 20 -np 5 [-c|-t] &")
        print("For qsub")
        print(f"{os.path.basename(__file__)} -e 'mllorenz.py' -js qsub -nch 12 -hl 7 -nn 15 -np 4 [-c|-t]] ")
        print("-c : check string, -t: short run for test algorithm")
        sys.exit(0)
    calc_ga(args.job_submit,args.exe,args.ifile,args.ofile,args.ndata,args.max_iter,args.ini_setup,args.nchromo,args.hlayers,args.nnodes,args.ngenerations,args.nparents_sets,args.prob_mutation, args.check_string, args.Ltest)
    return 0


if __name__ == '__main__':
    main()

