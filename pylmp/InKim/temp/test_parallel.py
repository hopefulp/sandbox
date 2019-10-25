import nutils as nu
import concurrent.futures
from functools import partial

def is_prime(n):
    if n < 2:
        return False
    if n is 2 or n is 3:
        return True
    if n % 2 is 0 or n % 3 is 0:
        return False
    if n < 9:
        return True
    k, l = 5, n ** 0.5
    while k <= l:
        if n % k is 0 or n % (k+2) is 0:
            return False
        k += 6
    return True

def process(n, r=10000):
    print("processing: {} ..< {}".format(n, n+r))
    s = sum((x for x in range(n, n+r) if is_prime(x) if x <= 2000000))
    #print(s)
    return s

@nu.timer
def main():
    r = 100000 
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as exe:
        result = 0
        for i in exe.map(partial(process, r=r), range(0, 2000000, r)):
            result += i
            #print(result)
        print(result)

main()
