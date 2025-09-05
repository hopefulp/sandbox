import numpy as np

def flat2to1( arr2d ):
    arr1d = arr2d.flatten()
    return arr1d

def extract_numeric_block(flines, i=0, j=-1):
    return [list(map(float, line.strip().split())) for line in flines[i:j]]

