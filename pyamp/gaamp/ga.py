import numpy
import random

# Selecting the best individuals in the current generation as parents for producing the offspring of the next generation.
selection_rule = 'best'

def select_mating_pool(chromo_mat, fit, num_parents):
    parents = numpy.empty((num_parents, chromo_mat.shape[1]))
    fitness = numpy.ndarray.tolist(fit)                     # change into list
    fitness.sort(reverse=True)                              # use decreasing order in case, sort()
    for i in range(num_parents):
        fitness_idx = numpy.where(fit == fitness[i])        # np.where returns tuple of array ( array([4,2,...), ... )
        parents[i] = chromo_mat[fitness_idx]
        #parents.append(chromo_mat[fitness_idx[0][0]])      # in case list
    return parents                                          

def crossover(parents, offspring_size):
    offspring = numpy.empty(offspring_size)
    # The point at which crossover takes place between two parents. Usually, it is at the center.
    crossover_point = numpy.uint8(offspring_size[1]/2)

    for k in range(offspring_size[0]):
        # Index of the first parent to mate.
        parent1_idx = k%parents.shape[0]
        # Index of the second parent to mate.
        parent2_idx = (k+1)%parents.shape[0]
        # The new offspring will have its first half of its genes taken from the first parent.
        offspring[k][0:crossover_point] = parents[parent1_idx][0:crossover_point]
        # The new offspring will have its second half of its genes taken from the second parent.
        offspring[k][crossover_point:] = parents[parent2_idx][crossover_point:]
    return offspring

def mutation(offspring_crossover, mutation_percent, nnode):
    num_mutations = numpy.uint8((mutation_percent*offspring_crossover.shape[1])/100)
    if num_mutations <= 0:
        num_mutations=1
    mutation_indices = numpy.array(random.sample(range(0, offspring_crossover.shape[1]), num_mutations))
    # Mutation changes a single gene in each offspring randomly.
    for idx in range(offspring_crossover.shape[0]):
        # The random value to be added to the gene.
        random_value = numpy.random.uniform(0, nnode, 1)//1
        if offspring_crossover[idx, mutation_indices] + random_value >nnode:
            offspring_crossover[idx, mutation_indices] = offspring_crossover[idx, mutation_indices] + random_value - nnode
        else:
            offspring_crossover[idx, mutation_indices] = offspring_crossover[idx, mutation_indices] + random_value  
    return offspring_crossover
