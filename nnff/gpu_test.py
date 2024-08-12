#!/home/joonho/anaconda3/bin/python

from torch import cuda

if cuda.is_available():
    print("cuda is available")


