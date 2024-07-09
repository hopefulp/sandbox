#!/home/joonho/anaconda3/bin/python
import numpy as np
import random
from common import whereami
Lprint = 0
# Selecting the best individuals in the current generation as parents for producing the offspring of the next generation.
selection_rule = 'best'


def select_mating_pool(chromo_mat, fitness, num_parents):
    '''
    obtaining new parents pool
        by score and num_parents
    '''
    print(f"in {whereami()} \n {chromo_mat} \n with fitness \n {fitness}")
    ### selection rule is in the order of best
    fsort = fitness.argsort()
    rank_descend = fsort[::-1][:num_parents]
    sel_fitness = fitness[rank_descend]
    sel_parents = chromo_mat[rank_descend]
    
    return sel_parents, sel_fitness


def crossover(parents, offspring_size):
    '''
    from parents make offspring using crossover
        offspring can be duplicated so all the offspring goes to mutate
    parents.shape = (nparents, nhl)
    offsprint_size = (nchromosome-nparents, nhl)
    '''
    offspring = np.empty(offspring_size, dtype=int)
    # The point at which crossover takes place between two parents. Usually, it is at the center.
    crossover_point = np.uint8(offspring_size[1]/2)
    check_parents_pairing = []
    for k in range(offspring_size[0]):
        ### get two parents, this might duplicate pairs
        ntry = 0
        while True:
            i = random.randint(0, parents.shape[0]-1)
            j = random.randint(0, parents.shape[0]-1)
            a = parents[i][:crossover_point]
            b = parents[j][crossover_point:]
            switch = random.randint(0, 1)
            if switch == 0:
                offspring[k] = np.concatenate((a,b))
            else:
                offspring[k] = np.concatenate((b,a))

            if [i, j, switch] in check_parents_pairing:
                ntry += 1
                if ntry > 5:
                    break
                continue
            else:
                check_parents_pairing.append([i, j, switch])
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
    if Lprint: print(f"num of mutations = {num_mutations} in {whereami()}")
    # Mutation changes only a single gene in each chromosome randomly.
    for ichromo in range(offspring_crossover.shape[0]):
        # The random value to be added to the gene.
        igene_mut = np.array(random.sample(range(0, offspring_crossover.shape[1]), num_mutations))
        nadd_gene = random.randint(0, nnode) # random.randint:[ ] -> 0 ~ nnode
        offspring_crossover[ichromo, igene_mut] = (offspring_crossover[ichromo, igene_mut] + nadd_gene) % (nnode+1) # to have 0~nnode

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
