'''
    Genetic Algorithm Package
'''
import sys
import numpy as np
import os
from .ga import *
import random
import time
import subprocess
import amp_ini
import pickle
from common import whereami, search_dirs

Lprint = 1

class GaAmp:
    '''
    checkfile_amptest  for checking test amp running in failed directory of 'amp-untrained-paramaters.amp'
    checkfile_generation to go to next generation when one generation is done, ampout_onegeneration_fitness='ga_gen_fit.txt'
    '''
    _chkfile_generation=amp_ini.ampout_onegeneration_fitness
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

    
    def run_generation(self, job_submit='qsub', dir_prefix='ch', istart=0, ngenerations=10, amp_job_setting=None, nparent_sets=4,
                        mutation_percent=10, check_print=False, sleeptime=10 ):
        cwd = os.getcwd()
        print(f"cwd {cwd} in {whereami()}")
        for igen in range(istart, ngenerations):
            print(f"{igen}-th generation in {ngenerations} in {whereami()}")
            #population = calc.chromo_mat.tolist() # Copy chromo_mat to population
            dir_genname = dir_prefix + '{:02d}'.format(igen) 
            fit = []                            # end of 1 generation
            gen_dirs=[]
            for i, chromo in enumerate(self.chromo_mat):
                chromo = [ x for x in chromo if x != 0] # remove 0 or HL
                if chromo == []:
                    chromo.append(1)
                dirname = dir_genname + '_{:02d}'.format(i)
                gen_dirs.append(dirname)
                os.system(f"make_dir.py {dirname} -w amp -j tr")

                ### submit nchromo-jobs inside each directory
                os.chdir(dirname)
                strHL=' '.join(str(x) for x in chromo)

                qjobname='ch{:02d}_{:02d}tr'.format(igen, i)
                if job_submit == 'qsub':
                    ampstr = amp_ini.Amp_string(queuename=qjobname, hl=strHL, add_amp_kw=amp_job_setting)
                    with open("mlet_tr.csh", 'w') as f:
                        f.write(ampstr.qscript)
                    print(ampstr())
                    os.system(ampstr())
                elif job_submit == 'node':
                    ampstr = amp_ini.Amp_string(hl=strHL, add_amp_kw=amp_job_setting)
                    str_node = "amp_run.py " + ampstr() + " -g &"
                    print(f"{str_node} in {whereami()}")
                    if check_print:
                        sys.exit(3)
                    os.system(str_node)

                os.chdir(cwd)
            ### to not use ''' ''' in the middle of function
            ### Wait On: until Training Step is done
            while True:
                time.sleep(sleeptime)
                ### Testing code: succeeded job, failed training jobs, and waiting
                dir_succ = search_dirs(dir_genname, 'amp.amp')
                dir_fail = search_dirs(dir_genname, 'amp-untrained-parameters.amp')
                ntraining = len(dir_succ) + len(dir_fail)
                nwait = self.nchromo - ntraining
                if nwait:
                    print(f"succeed dir {len(dir_succ)} failed dir {len(dir_fail)}, waiting = {nwait} for {sleeptime} sec")
                else:
                    break
            ### TEST direct run for test for all directories
            for test_dir in gen_dirs:  # run test in all dir
                os.chdir(test_dir)
                dir_suf=test_dir[-2:]
                qjobname='ch{:02d}_{:2}te'.format(igen,dir_suf)
                if job_submit == 'qsub':
                    ### set test here
                    amp_job_setting['ampjob'] = 'gate'
                    amp_job_setting['ncore'] = '1'
                    amp_job_setting['max_iter'] = None
                    amp_job_setting['elimit']= None
                    amp_job_setting['mem'] = '3G'
                    ampstr = amp_ini.Amp_string(quenename=qjobname, hl=strHL, add_amp_kw=amp_job_setting)
                    with open("mlet_te.csh", 'w') as f:
                        f.write(ampstr.qscript)
                    print(ampstr())
                elif job_submit == 'node':
                    ### run test force
                    ampstr = amp_ini.Amp_string(jsubmit=amp_job_setting['jsubmit'], ampjob='tega', dlist=amp_job_setting['dlist'])
                    str_node = "amp_run.py " + ampstr() + " -g &"
                    print(str_node)
                    os.system(str_node)
                os.chdir(cwd)
            ### Wait On: until one generation is done defined by number of lines in GaAmp._chkfile_generation
            while True:
                time.sleep(sleeptime)
                pftxt = os.path.isfile(GaAmp._chkfile_generation)
                if pftxt == True:
                    nline = subprocess.check_output(['wc','-l',GaAmp._chkfile_generation])
                    n = eval(nline.decode('utf-8').split()[0])
                    print(f"n line in {GaAmp._chkfile_generation} = {n}")
                    if n == self.nchromo:
                        with open(GaAmp._chkfile_generation, "r") as f:
                            ga_gen_fit = f.read().splitlines()  # remove \n character
                        os.system(f"rm {GaAmp._chkfile_generation}") # this file transfered data to alines so removed for next generation
                        break
                else:
                    print(f'{GaAmp._chkfile_generation} is not made yet')
                
            ### analyze the input of ga_gen_fit.txt
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
           
            with open("GANNtest"+cwd[-2:]+"_Genes.txt", "a") as f:
                f.write(f"Generation : {igen}\nGenes : \n{gen_chromos}\nfitness : \n{gen_fitnesses}\n")
            with open("GANNtest"+cwd[-2:]+".txt", "a") as f:
                f.write(f"Generation : {igen}\nGenes : {gen_chromos[i_bestfit]}, fitness : {best_fit}\n")
                
            self.chromo_mat = gen_step(gen_chromos, gen_fitnesses, nparent_sets, best_fit, mutation_percent, self.nnode)
            print(f"go to next generation {self.chromo_mat}")
            fit=[]

def gen_step(chromo_mat, fit, nparent_sets, best_fit, mutation_percent, nnode):
    parents = select_mating_pool(chromo_mat, fit, nparent_sets)     # returns list
    print(f"parents = {parents}")   

    offspring_crossover = crossover(parents, offspring_size=(chromo_mat.shape[0]-parents.shape[0], chromo_mat.shape[1]))
    print(f"offspring_crossover = {offspring_crossover}")
    offspring_mutation = mutation(offspring_crossover, mutation_percent, nnode)
    print(f"offspring_mutation = {offspring_mutation}")
    chromo_mat[0:nparent_sets] = parents
    chromo_mat[nparent_sets:] = offspring_crossover

    return chromo_mat
