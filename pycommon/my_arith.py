from math import factorial
import numpy as np

### product all the elements in 1D array
def product(iterable):
    prod = 1
    for n in iterable:
        prod *= n
    return prod

def npr(n, r):
    """
    Calculate the number of ordered permutations of r items taken from a
    population of size n.

    >>> npr(3, 2)
    6
    >>> npr(100, 20)
    1303995018204712451095685346159820800000
    """
    assert 0 <= r <= n
    return product(range(n - r + 1, n + 1))

def ncr(n, r):
    """
    Calculate the number of unordered combinations of r items taken from a
    population of size n.

    >>> ncr(3, 2)
    3
    >>> ncr(100, 20)
    535983370403809682970
    >>> ncr(100000, 1000) == ncr(100000, 99000)
    True
    """
    assert 0 <= r <= n
    if r > n // 2:
        r = n - r
    return npr(n, r) // factorial(r)

def divide_int(inum, n):
    ave = inum / float(n)
    li = [ int(i*ave) for i in range(n) if i != 0]
    li.append(inum-1)
    return li

def rmse(y, y_bar):
    err = np.sqrt(np.square(np.array(y) - np.array(y_bar)).mean())
    return err
