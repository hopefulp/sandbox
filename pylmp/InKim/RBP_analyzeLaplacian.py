import numpy as np
import networkx as nx

def find_smallest_eigval(L):
    """
    find the second smallest eigenvalue (the smallest one except 0) from the list of eigvals of Laplacian L
    """
    l = [i for i in L if i > 1e-10]
    return min(l)
