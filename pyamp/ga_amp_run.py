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

def calc_ga(fin, job_submit, setup, nchromo, nhl, nnode, ngenerations, nparent_sets, mutation_percent):

    calc = GaAmp(setup=setup, nchromo=nchromo, nhl=nhl, nnode=nnode)

    ### choose whether train force [any True] or None
    force_train = None      # for force_train 1, for energy only None
    sleep_time  = 10        # 10 for test, 300 for run

    start_point = 0
    subdir_prefix = 'ch'

    #p = Parameters({'jsubmit':job_submit,'ampjob':'trga', 'elimit':'0.0001', 'ncore':'8', 'max_iter':500,
    #                'mem':'12G', 'dlist':[1000,2000,3500,3600]})
    ### for Test
    p = Parameters({'jsubmit':job_submit,'ampjob':'trga', 'elimit':'0.0001', 'ncore':'2', 'max_iter':100,
                    'mem':'12G', 'dlist':[1000,1100,3500,3510]})
    ### include force or not
    if force_train == None:     # else default p.train_f is set
        p.train_f = None

    # make jobstring and pass to calc or pass parameters then calc.run_generation make the jobstring due to hl made in run_generation
    check_ampstring = False  # do not run job

    calc.run_generation(job_submit=job_submit, dir_prefix=subdir_prefix, istart=start_point, ngenerations=ngenerations, 
                        amp_job_setting=p, nparent_sets=nparent_sets, 
                        mutation_percent=mutation_percent, check_print=check_ampstring, sleeptime=10)

    print('GA CALCULATION IS DONE')
    return 0

def main():
    parser = argparse.ArgumentParser(description='Process for Genetic Algorithm & Artificial Neural Network')
    parser.add_argument('-f', '--inf', default='OUTCAR', help='input file')
    parser.add_argument('-js', '--job_submit', default='qsub', choices=['qsub','node'], help='Select between node and login server')
    parser.add_argument('-i', '--ini_setup', default = 'random', choices=['random','reference','ongoing'], help='reference/ongoing from GANNtest file')
    parser.add_argument('-nch', '--nchromo', default=20, type=int, help="number of HLs in each generation")
    parser.add_argument('-hl', '--hlayers', default=10, type=int, help="fixed number of layers in each HL")
    parser.add_argument('-nn', '--nnodes', default=15, type=int, help="number of nodes in each HL")
    parser.add_argument('-np', '--nparents_sets', default=4, type=int, help="number of parents mating: select best hls")
    parser.add_argument('-mp', '--prob_mutation', default=10, type=int, help="percentage of mutation")
    parser.add_argument('-ng', '--ngenerations', default=10, type=int, help="number of generation")
    args = parser.parse_args()

    pf = glob.glob("ga_fit.txt")
    if pf != []:
        os.system("rm ga_fit.txt")
    #fin = glob.glob("*extxyz") 
    calc_ga(args.inf,args.job_submit,args.ini_setup,args.nchromo,args.hlayers,args.nnodes,args.ngenerations,args.nparents_sets,args.prob_mutation)
    return 0


if __name__ == '__main__':
    main()

