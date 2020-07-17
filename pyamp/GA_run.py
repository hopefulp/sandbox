#!/home/joonho/anaconda3/bin/python

import numpy as np
import os
import GA
import random
import shutil
import glob
import time
import argparse
import qsub_restart

def ran_pop(num_pop, layers, nodes, genes, initial_pop, cnt):
    for curr_sol in np.arange(0, num_pop): # 30 Chromosomes have 5 genes
        for layer in layers:
            genes.append(random.choice(nodes))
        initial_pop.append(genes)
        genes=[]

    for pop in initial_pop:
        if 0 in pop:
            cnt = 0
        while 0 in pop:
            pop.remove(0)
            cnt += 1
        for i in range(cnt):
            pop.append(0)
            if i == cnt-1:
                cnt = 0
        if pop == [0]*len(layers):
            pop.append(1)

    return initial_pop

def desig_pop(cwd, num_pop, initial_pop, layers):
    len_layer = len(layers)
    with open('GANNtest'+cwd+'_Genes.txt','r') as f:
        flines = f.readlines()
    for i in range(len(flines)):
        if 'Generation' in flines[-(i+1)]:
            break
    generation = eval(flines[-(i+1)][-3:])
 
        
    pop_lines = []
    tmp_list = []
    pop_lines.append(flines[-(i+1)+2:-(i+1)+2+num_pop])
    for j in range(len(pop_lines)):
        for s in [',',"'",'[',']','\n','\\n']:
            pop_lines[j] = str(pop_lines[j])
            pop_lines[j] = pop_lines[j].replace(s,'')
        pop_lines[j] = pop_lines[j].split(' ')

        while '' in pop_lines[j]:
            pop_lines[j].remove('')
    
        cnt = 0
        for k in pop_lines[j]:
            cnt += 1
            tmp_list.append(eval(k))
            if cnt == len_layer :
                initial_pop.append(tmp_list)
                tmp_list = []
                cnt = 0

    for i in range(len(flines)):
        if 'Error' in flines[-(i+1)]:
            break

    fit_lines = flines[-(i+1)+1:]
    fitness = []
    for j in range(len(fit_lines)):
        for s in ['[',']','\n']:
            fit_lines[j] = str(fit_lines[j])
            fit_lines[j] = fit_lines[j].replace(s,'')
        fit_lines[j] = fit_lines[j].split(' ')
       
        while '' in fit_lines[j]:
            fit_lines[j].remove('')

        for k in fit_lines[j]:
            tmp_list.append(eval(k))
        fitness.extend(tmp_list)
        tmp_list = []

    return (generation, initial_pop, fitness) 
        
def gen_step(pop_mat, fit, num_parents_mating, best_fit, mutation_percent, max_node, min_node): 
    parents = GA.select_mating_pool(pop_mat, fit, num_parents_mating)

    offspring_crossover = GA.crossover(parents, offspring_size=(pop_mat.shape[0]-parents.shape[0], pop_mat.shape[1]))

    offspring_mutation = GA.mutation(offspring_crossover, mutation_percent, max_node, min_node)

    pop_mat[0:num_parents_mating] = parents
    pop_mat[num_parents_mating:] = offspring_crossover
 
    return pop_mat


def calc_ga(fin, init_pop, num_pop, max_nhl, nnodes, num_parents_mating, mutation_percent, num_generations, queue_job='qsub'):
    cwd=os.getcwd()
    cwd = cwd[-2:] 
   
    #num_pop = 20                        
    #num_parents_mating = 4         # select best hls
    #num_generations = 10
    #mutation_percent = 10
    layers = np.arange(max_nhl)  
    nodes = np.arange(nnodes)
    max_node = max(nodes)
    min_node = min(nodes)
    initial_pop=[]
    genes = []
    fit = []
    pop_mat = []   
 
    start_point = 0
    cnt = 0
        
    if init_pop == 'random': 
        pop_mat = ran_pop(num_pop, layers, nodes, genes, initial_pop, cnt)
        pop_mat = np.array(pop_mat)
        for i in range(num_pop):
           dir_name='./pop{:02d}'.format(i)
           copy_name='./pop{:02d}/.'.format(i)
           os.system("rm -rf "+dir_name)
           os.mkdir(dir_name)
           #os.system("cp "+fin[0]+" amp_ga.py my* GA.py "+copy_name)    # my* ?
           os.system(f"ln -s {fin[0]} {copy_name}/{fin[0]}" 

    elif init_pop == 'reference':
        option = desig_pop(cwd, num_pop, initial_pop, layers)
        pop_mat = option[1]
        pop_mat = np.array(pop_mat)
        start_point = option[0]+1
        fitness = option[2]
        fitness = np.array(fitness)
        best_fit = min(fitness)
        pop_mat = gen_step(pop_mat, fitness, num_parents_mating, best_fit, mutation_percent, max_node, min_node)
        np.array(pop_mat)
    
    elif init_pop == 'ongoing':
        while True:
            time.sleep(10)
            pftxt = os.path.isfile('pop_fit.txt')
            if pftxt == True:
                alines=[]
                with open("pop_fit.txt", "r") as a:
                    alines=a.readlines()
                    if len(alines) >= num_pop:
                        print('delete pop_fit')
                        os.system("rm pop_fit.txt")
                        break
                    elif len(alines) < num_pop:
                        if queue_job == 'qsub':
                            qsub_restart.qsub(num_pop)

        pop_tmp = []                                           
        for i in range(len(alines)):
            tmp_list = alines[i].split(']')
            pop_list = []
            
            for s in ['[',' ']:
                tmp_list[0] = tmp_list[0].replace(s,'')

            tmp_list[0] = tmp_list[0].split(',')             
            
            for l in tmp_list[0]:
                pop_list.append(eval(l))

            if len(pop_list) != len(layers):
                for cnt in range(len(layers)-len(pop_list)):
                    pop_list.append(0)         
            pop_tmp.append(pop_list)
            fit.append(eval(tmp_list[1]))

        pop_mat = np.array(pop_tmp)
        fit = np.array(fit)
        best_fit = min(fit)
 
        with open('GANNtest'+cwd+'_Genes.txt','r') as f:
            flines = f.readlines()
        for i in range(len(flines)):
            if 'Generation' in flines[-(i+1)]:
                break
        generation = eval(flines[-(i+1)][-3:])

      
        with open("GANNtest"+cwd+"_Genes.txt", "a") as f:
            f.write("Generation : {}\nGenes : \n{}\nError : \n{}\n".format(generation, pop_mat, fit))
        with open("GANNtest"+cwd+".txt", "a") as f:
            f.write("Generation : {}\nGenes : {}, Error : {}\n".format(generation, pop_mat[np.where(fit == best_fit)], best_fit))
            
        pop_mat = gen_step(pop_mat, fit, num_parents_mating, best_fit, mutation_percent, max_node, min_node)
        fit=[]
 
        option = desig_pop(cwd, num_pop, initial_pop, layers)
        pop_mat = option[1]
        pop_mat = np.array(pop_mat)
        start_point = option[0]+2
        fitness = option[2]
        fitness = np.array(fitness)
        best_fit = min(fitness)
        pop_mat = gen_step(pop_mat, fitness, num_parents_mating, best_fit, mutation_percent, max_node, min_node)
        np.array(pop_mat)
 
    for generation in range(start_point, num_generations):
        population = np.ndarray.tolist(pop_mat) # Copy pop_mat to population
        cnt = 0
        for pop in population:
            while 0 in pop:
                pop.remove(0) # Because amp cannot recognize 0 node
            
            if pop == []:
                pop.append(1)
            with open('./pop{:02d}/HL'.format(cnt),'w') as f:
                f.write('{}'.format(pop)) 
            strHL=' '.join(str(x) for x in pop)
  
            with open('./pop{:02d}/amp_pot.sh'.format(cnt), 'w') as f:
                f.write('#!/usr/bin/csh\n\n')
                f.write('$HOME/anaconda3/bin/python amp_ga.py '+fin[0]+' train -n 750 -hl '+strHL+' -el 0.0001\n')
            
            os.chdir('./pop{:02d}'.format(cnt))
            if os.path.isfile('score') == True:
                os.remove('score')
            rmfin = glob.glob('amp_pot.sh.e*')
            for rm in rmfin:            
                os.remove(rm)
            
            if queue_job == 'qsub':
                os.system('qsub -cwd amp_pot.sh') # Need to change this command. csh -> qsub -cwd 
            elif queue_job == 'csh':
                os.system('csh amp_pot.sh &')
                os.system("pidof amp_pot.sh > pid")
                with open('pid','r') as f:
                    flines = f.readlines()
                flines = ' '.join(flines)
                flines = flines.split(' ')
                flines[-1]
                with open('pid','w') as f:
                    f.write(flines[-1])

            os.chdir('../')
            cnt += 1
            
            # From this line, code how to extract information from pop_fit.txt
           
        while True:
            time.sleep(10)
            pftxt = os.path.isfile('pop_fit.txt')
            if pftxt == True:
                alines=[]
                with open("pop_fit.txt", "r") as a:
                    alines=a.readlines()
                    if len(alines) >= num_pop:
                        os.system("rm pop_fit.txt")
                        break
                    elif len(alines) < num_pop:
                        if queue_job == 'qsub':
                            qsub_restart.qsub(num_pop)
                        elif queue_job == 'csh':
                            qsub_restart.csh(num_pop)
        pop_tmp = []
        for i in range(len(alines)):
            tmp_list = alines[i].split(']')
            pop_list = []
            
            for s in ['[',' ']:
                tmp_list[0] = tmp_list[0].replace(s,'')

            tmp_list[0] = tmp_list[0].split(',')             
            
            for l in tmp_list[0]:
                pop_list.append(eval(l))

            if len(pop_list) != len(layers):
                for cnt in range(len(layers)-len(pop_list)):
                    pop_list.append(0)         
            pop_tmp.append(pop_list)
            fit.append(eval(tmp_list[1]))

        pop_mat = np.array(pop_tmp)
        fit = np.array(fit)
        best_fit = min(fit)
       
        with open("GANNtest"+cwd+"_Genes.txt", "a") as f:
            f.write("Generation : {}\nGenes : \n{}\nError : \n{}\n".format(generation, pop_mat, fit))
        with open("GANNtest"+cwd+".txt", "a") as f:
            f.write("Generation : {}\nGenes : {}, Error : {}\n".format(generation, pop_mat[np.where(fit == best_fit)], best_fit))
            
        pop_mat = gen_step(pop_mat, fit, num_parents_mating, best_fit, mutation_percent, max_node, min_node)
        fit=[]

    '''
        parents = GA.select_mating_pool(pop_mat, fit, num_parents_mating)

        offspring_crossover = GA.crossover(parents, offspring_size=(pop_mat.shape[0]-parents.shape[0], pop_mat.shape[1]))

        offspring_mutation = GA.mutation(offspring_crossover, mutation_percent, max_node, min_node)
  
        pop_mat[0:num_parents_mating] = parents
        pop_mat[num_parents_mating:] = offspring_crossover
    '''
        
    return print('ANN CALCULATION IS DONE')

def main():
    parser = argparse.ArgumentParser(description='Process for Genetic Algorithm & Artificial Neural Network')
    parser.add_argument('-i', '--init_pop', default = 'random', choices=['random','reference','ongoing'], help='reference/ongoing from GANNtest file')
    #parser.add_argument('-qj', '--queue_job', default='qsub', choices=['qsub','csh'], help='Select between node and login server')
    parser.add_argument('-f', '--inf', default=['OUTCAR'], nargs='+', help='input file')
    parser.add_argument('-nhls', '--hl_nsets', default=20, type=int, help="number of HLs in each generation")
    parser.add_argument('-nhl', '--nhlayers', default=10, type=int, help="fixed number of layers in each HL")
    parser.add_argument('-nn', '--nnodes', default=15, type=int, help="number of nodes in each HL")
    parser.add_argument('-npm', '--nparents_mate', default=4, type=int, help="number of parents mating")
    parser.add_argument('-mp', '--pmutation', default=10, type=int, help="percentage of mutation")
    parser.add_argument('-ng', '--ngenerations', default=5, type=int, help="number of generation")
    args = parser.parse_args()

    pf = glob.glob("pop_fit.txt")
    if pf != []:
        os.system("rm pop_fit.txt")
    #fin = glob.glob("*extxyz") 
    calc_ga(args.inf, args.init_pop, args.hl_nsets, args.nhlayers, args.nnodes, args.nparents_mate, args.pmutation, args.ngenerations)
    return 0


if __name__ == '__main__':
    main()

