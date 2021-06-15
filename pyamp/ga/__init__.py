'''
    Genetic Algorithm Package
'''
import sys
import numpy as np
import os
import random
import time
import ga_ini
import subprocess
import pickle
from common import whereami, search_dirs

from . import ga_sub as ga

Lprint = 1

class GA:
    '''
    checkfile_amptest
    checkfile_generation to go to next generation when one generation is done, fitness_onegen='ga_gen_fit.txt'
    '''
    _chkfile_generation=ga_ini.fitness_onegen
    def __init__(self, setup='random', nchromo=10, nhl=10, nnode=15):
        self.setup      = setup
        self.nchromo    = nchromo
        self.nhl        = nhl
        self.nnode      = nnode
        if self.setup == 'random':
            self.chromo_mat = self.random_population()
        else:
            print(f"Error:: can't set self.chromo_mat {whereami()}")
            sys.exit(10)


    def random_population(self):
        '''
        randomly generate population matrix
        int nchromo: number of population of chromosome as row in matrix
        int nhl: number of max layers
        '''
        #print(f"nlyaers = {nlayers} in {whereami()}")
        chromo_mat=[]
        for curr_sol in np.arange(0, self.nchromo): # 
            chromosome=[]
            for i in range(self.nhl):
                chromosome.append(random.randint(0, self.nnode))
            chromo_mat.append(chromosome)
        if Lprint: print(f"chromosom_matrix = \n{np.array(chromo_mat)}\n in {whereami()}")
        ### move 0 to the end and make nlayers chromo
        for chromo in chromo_mat:
            while 0 in chromo:
                chromo.remove(0)
            for i in range(self.nhl-len(chromo)):
                chromo.append(0)
            if chromo == [0]*self.nhl:
                chromo[0]=1             # not to make [0,0...]
        return chromo_mat

    def design_pop(cwd, nchromo, initial_pop, layers):
        '''
        in case ini_setup == 
        designn population matrix
        '''
        len_layer = len(layers)
        with open('GANNtest'+cwd+'_Genes.txt','r') as f:
            flines = f.readlines()
        for i in range(len(flines)):
            if 'Generation' in flines[-(i+1)]:
                break
        generation = eval(flines[-(i+1)][-3:])
     
            
        pop_lines = []
        tmp_list = []
        pop_lines.append(flines[-(i+1)+2:-(i+1)+2+nchromo])
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

    
    def run_generation(self, job_submit='qsub', fname=None, dir_prefix='ch', istart=0, ngenerations=10, ml_job_setting=None,\
                        nparent_sets=4, mutation_percent=10, check_print=False, sleeptime=10 ):
        cwd = os.getcwd()
        print(f"job_submit {job_submit} in {whereami()}")
        chromo_parents=[]
        for igen in range(istart, ngenerations):
            print(f"{igen}-th generation in {ngenerations} in {whereami()}")
            ### need to remove GA._chkfile_generation used for indication of cal_generation finished
            if os.path.isfile(GA._chkfile_generation):
                os.system(f"rm {GA._chkfile_generation}")
            dir_genname = dir_prefix + '{:02d}'.format(igen) 
            fit = []                            # end of 1 generation
            gen_dirs=[]
            for i, chromo in enumerate(self.chromo_mat):
                chromo = [ int(x) for x in chromo if x != 0]         # chromo is contracted
                if chromo == []:
                    chromo.append(1)
                dirname = dir_genname + '_{:02d}'.format(i)
                gen_dirs.append(dirname)
                ###### in case make dir
                #os.system(f"make_dir.py {dirname} -w amp -j tr")
                ### submit nchromo-jobs inside each directory
                #os.chdir(dirname)
                strHL=' '.join(str(x) for x in chromo)

                qjobname='ch{:02d}_{:02d}tr'.format(igen, i)
                if job_submit == 'qsub':
                    gastr = ga_ini.Ga_string(queuename=qjobname,jsubmit=job_submit, hl=strHL, add_ml_kw=ml_job_setting)
                    with open("mlet_tr.csh", 'w') as f:
                        f.write(gastr.qscript)
                    str_2run = gastr()
                ### node includes 'gpu'
                elif job_submit == 'node':
                    gastr = ga_ini.Ga_string(jsubmit=job_submit, hl=strHL, add_ml_kw=ml_job_setting)
                    str_2run = gastr() + f' -ig {igen} -ga'
                print(f"{str_2run} in {whereami()}")
                if check_print:
                    sys.exit(3)
                os.system(str_2run)

                #os.chdir(cwd)
            ### to not use ''' ''' in the middle of function
            ### Wait On: until Training Step is done
            ### One generation is done in "while"
            while True:
                time.sleep(sleeptime)
                ### Testing code: succeeded job, failed training jobs, and waiting
                results = subprocess.check_output(f'wc -l {GA._chkfile_generation}', shell=True) 
                wl = int(results.split()[0])
                nwait = self.nchromo - wl
                if nwait:
                    print(f"done {wl}, waiting = {nwait} for {sleeptime} sec in {igen}-gen")
                else:
                    break
            ### TEST direct run for test for all directories
            ### summary one generatioin
    ### for igen in range(istart, ngenerations):
            ### analyze ga_gen_fit.txt of GA-one generation output
            print("one generation is done")
            with open(GA._chkfile_generation, "r") as f:
                ga_gen_fit = f.read().splitlines()
            gen_chromos=[]
            gen_fitnesses=[]
            for hl_fit in ga_gen_fit:
                print(hl_fit)
                hllist = hl_fit.split()
                fitness = eval(hllist.pop())    # extract last element then evaluate string to float
                hl = list(map(int, hllist))     # map change each elements then change obj to int
                hl.extend([0] * (self.nhl-len(hllist))) # hl reserve the original size of chromo_matrix 
                
                gen_chromos.append(hl)
                gen_fitnesses.append(fitness)

            gen_chromos     = np.array(gen_chromos)
            print(f"result  = {gen_chromos.shape}")
            gen_fitnesses   = np.array(gen_fitnesses)
            best_fit        = np.max(gen_fitnesses)
            i_bestfit       = np.argmax(gen_fitnesses)

            ### append parents to gen_chromos to select parents for next generation and fitness also
            if chromo_parents == []:
                chromo_parents = np.zeros((nparent_sets, gen_chromos.shape[1]))
                fit_parents    = np.zeros((nparent_sets))
            print(f"size of {chromo_parents}, {gen_chromos}")
            gen_chromos_par     = np.vstack((chromo_parents, gen_chromos))
            gen_fitnesses_par   = np.append(fit_parents, gen_fitnesses)
            print(f"{gen_chromos_par} {gen_fitnesses_par}")

            with open("GANNtest"+cwd[-2:]+"_all.txt", "a") as f:
                f.write(f"Generation : {igen}\nGenes : \n{gen_chromos_par}\nfitness : \n{gen_fitnesses_par}\n")
            with open("GANNtest"+cwd[-2:]+"_best.txt", "a") as f:
                f.write(f"Generation : {igen}\nGenes : {gen_chromos[i_bestfit]}, fitness : {best_fit}\n")

            ### reduce dimension of gen_chromos if necessary
            gen_chromos_par = gen_chromos_par[:, ~np.all(gen_chromos==0, axis=0)]   # reduce hl-layers if col==0
            if self.nhl != gen_chromos_par.shape[1]:
                print("reduce size of rank-1")
                self.nhl = gen_chromos_par.shape[1]

            self.chromo_mat, chromo_parents, fit_parents = gen_step(gen_chromos_par, gen_fitnesses_par, nparent_sets, best_fit, mutation_percent, self.nnode)
            print(f"{chromo_parents} .. {fit_parents}")
            print(f"go to next generation {self.chromo_mat}")
            fit=[]


def gen_step(chromo_mat, fit, nparent_sets, best_fit, mutation_percent, nnode):
    '''
    chromo_mat includes nparents which will not run
    chromo_mat.size = self.chromo_mat.size(offspring==running jobs) + nparents
    '''
    parents, fit_parents = ga.select_mating_pool(chromo_mat, fit, nparent_sets)     # returns list
    print(f"parents = {parents}")   

    offspring_crossover = ga.crossover(parents, offspring_size=(chromo_mat.shape[0]-nparent_sets, chromo_mat.shape[1]))
    print(f"offspring_crossover = {offspring_crossover}")
    offspring_mutation = ga.mutation(offspring_crossover, mutation_percent, nnode)
    print(f"offspring_mutation = {offspring_mutation}")
    #chromo_mat[0:nparent_sets] = parents
    #chromo_mat[nparent_sets:] = offspring_crossover

    return offspring_mutation, parents, fit_parents

