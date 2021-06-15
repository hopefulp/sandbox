#!/home/joonho/anaconda3/bin/python
import numpy as np
import random
from common import whereami
Lprint = 1
# Selecting the best individuals in the current generation as parents for producing the offspring of the next generation.
selection_rule = 'best'

def select_mating_pool(chromo_mat, fitness, num_parents):
    print(f"{chromo_mat} {fitness} in {whereami()}")
    chromo_size = chromo_mat.shape[1]
    #parents = np.empty((num_parents, chromo_size))
    fitness_list = fitness.tolist()                                 # change into list
    fitness_list.sort(reverse=True)                                 # list is sorted
    sel_fitness = fitness_list[:num_parents]
    par=[]
    for i in range(num_parents):
        fitness_idx = np.where(fitness == fitness_list[i])          # np.where returns tuple of array ( array([4,2,...), ... )
        par.append(chromo_mat[fitness_idx].reshape((chromo_size,)))         # parents[i] = chromo_mat[idx] makes int to float why?
    print(f"{np.array(par)} {sel_fitness} in {whereami()}")        
    return np.array(par), np.array(sel_fitness)                                       # return np.array


def crossover(parents, offspring_size):
    '''
    from parents make offspring using crossover
        offspring can be duplicated so all the offspring goes to mutate
    parents.shape = (nparents, nhl)
    offsprint_size = (nchromosome-nparents, nhl)
    '''
    offspring = np.empty(offspring_size)
    # The point at which crossover takes place between two parents. Usually, it is at the center.
    crossover_point = np.uint8(offspring_size[1]/2)
    pair_parent = []
    for k in range(offspring_size[0]):
        ### get two parents, this might duplicate pairs
        ntry = 0
        while True:
            i = random.randint(0, parents.shape[0]-1)
            j = random.randint(0, parents.shape[0]-1)
            offspring[k][0:crossover_point] = parents[i][0:crossover_point]
            offspring[k][crossover_point:] = parents[j][crossover_point:]
            if [i, j] in pair_parent:
                ntry += 1
                if ntry > 5:
                    break
                continue
            else:
                pair_parent.append([i, j])
                break
    return offspring

def mutation(offspring_crossover, mutation_percent, nnode):
    '''
    offspring_crossover: all the offspring
        all the offspring will mutate
        do not run parents again
    '''
    num_mutations = np.uint8((mutation_percent*offspring_crossover.shape[1])/100)
    if num_mutations <= 0:
        num_mutations=1
    igene_mut = np.array(random.sample(range(0, offspring_crossover.shape[1]), num_mutations))
    # Mutation changes only a single gene in each chromosome randomly.
    for ichromo in range(offspring_crossover.shape[0]):
        # The random value to be added to the gene.
        nadd_gene = np.random.randint(0, nnode-1) # if nnode is added, it will not change
        if (offspring_crossover[ichromo, igene_mut] + nadd_gene) > nnode:
            offspring_crossover[ichromo, igene_mut] = offspring_crossover[ichromo, igene_mut] + nadd_gene - nnode
        else:
            offspring_crossover[ichromo, igene_mut] = offspring_crossover[ichromo, igene_mut] + nadd_gene  
    return offspring_crossover

def random_pop(nchromo, nhl, nnode):
    '''
    copied from gaamp.__init__
    randomly generate population matrix
    int nchromo: number of population of chromosome as row in matrix
    int nhl: number of max layers
    '''
    #print(f"nlyaers = {nlayers} in {whereami()}")
    chromo_mat=[]
    for curr_sol in np.arange(0, nchromo): # 
        chromosome=[]
        for i in range(nhl):
            chromosome.append(random.randint(0, nnode))
        chromo_mat.append(chromosome)
    if Lprint: print(f"chromosom_matrix = \n{np.array(chromo_mat)}\n in {whereami()}")
    ### move 0 to the end and make nlayers chromo
    for chromo in chromo_mat:
        while 0 in chromo:
            chromo.remove(0)
        for i in range(nhl-len(chromo)):
            chromo.append(0)
        if chromo == [0]*nhl:
            chromo[0]=1             # not to make [0,0...]
    return np.array(chromo_mat)



def main():
    ''' test functions'''
    nchromo=5
    chromo_mat = random_pop(nchromo, 5, 10)
    fitness = np.random.rand(nchromo) * 10
    parents, fitness = select_mating_pool(chromo_mat, fitness, 3)
    print(f"{parents} {fitness} in {whereami()}")
    return 0

if __name__ == '__main__':
    main()
