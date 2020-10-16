#!/home/joonho/anaconda3/bin/python

import numpy as np
import os
from gaamp import GaAmp
import glob
import argparse
#from gaamp import qsub_restart
import sys
from common import whereami 
from ase.calculators.calculator import Parameters

def calc_ga(fin, job_submit, setup, nchromo, nhl, nnode, ngenerations, nparent_sets, mutation_percent, Lcheck_string, Ltest):

    calc = GaAmp(setup=setup, nchromo=nchromo, nhl=nhl, nnode=nnode)
    ### whether train force or not
    force_train = 1         # [1,any|None], True: force_train ; None: energy-only train

    ### run_generation params
    start_point = 0
    subdir_prefix = 'ch'

    if Ltest:
        sleeptime = 10
        p = Parameters({'jsubmit':job_submit,'ampjob':'trga', 'elimit':'0.0001', 'ncore':'8', 'max_iter':20,
                    'mem':'5G', 'dlist':[1000,1100,3500,3510]})
    else:
        sleeptime = 300
        p = Parameters({'jsubmit':job_submit,'ampjob':'trga', 'elimit':'0.0001', 'ncore':'8', 'max_iter':500,
                    'mem':'12G', 'dlist':[1000,2000,3500,3600]})
    ### include force or not
    if force_train == None:     # else default p.train_f is set
        p.train_f = None

    # make jobstring and pass to calc or pass parameters then calc.run_generation make the jobstring due to hl made in run_generation
    calc.run_generation(job_submit=job_submit, dir_prefix=subdir_prefix, istart=start_point, ngenerations=ngenerations, 
                        amp_job_setting=p, nparent_sets=nparent_sets, 
                        mutation_percent=mutation_percent, check_print=Lcheck_string, sleeptime=sleeptime)

    print('GA CALCULATION IS DONE')
    return 0

def main():
    parser = argparse.ArgumentParser(description='Process for Genetic Algorithm & Artificial Neural Network')
    parser.add_argument('-f', '--inf', default='OUTCAR', help='input file')
    parser.add_argument('-js', '--job_submit', default='qsub', choices=['qsub','node'], help='Select between node and login server')
    parser.add_argument('-i', '--ini_setup', default = 'random', choices=['random','reference','ongoing'], help='reference/ongoing from GANNtest file')
    parser.add_argument('-u', '--usage', action='store_true', help='show usage')
    parser.add_argument('-c', '--check_string', action='store_true', help='check command string of amp_run.py but not submit the jobs')
    parser.add_argument('-t', '--Ltest', action='store_true', help='run test job for checking source code')
    genetic_params = parser.add_argument_group(title = 'Genetic Algorithm')
    genetic_params.add_argument('-nch', '--nchromo', default=20, type=int, help="number of HLs in each generation")
    genetic_params.add_argument('-hl', '--hlayers', default=10, type=int, help="fixed number of layers in each HL")
    genetic_params.add_argument('-nn', '--nnodes', default=15, type=int, help="number of nodes in each HL")
    genetic_params.add_argument('-np', '--nparents_sets', default=4, type=int, help="number of parents mating: select best hls")
    genetic_params.add_argument('-mp', '--prob_mutation', default=10, type=int, help="percentage of mutation")
    genetic_params.add_argument('-ng', '--ngenerations', default=10, type=int, help="number of generation")

    args = parser.parse_args()
    if args.usage:
        print("For Node")
        print("ga_amp_run.py -js node -nch 8 -hl 10 -nn 10 -np 5 [-c|-t] &")
        print("For qsub")
        print("ga_amp_run.py -js qsub -nch 12 -hl 7 -nn 15 -np 4 [-c|-t]] ")
        print("-c : check string, -t: short run for test algorithm")
        sys.exit(0)

    calc_ga(args.inf,args.job_submit,args.ini_setup,args.nchromo,args.hlayers,args.nnodes,args.ngenerations,args.nparents_sets,args.prob_mutation, args.check_string, args.Ltest)
    return 0


if __name__ == '__main__':
    main()

