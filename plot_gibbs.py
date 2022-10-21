#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os

def plot_gibbs(fname):
    ene=[]
    with open(fname, 'r') as f:
        lines=f.readlines()
        for line in lines:
            en = float(line.strip())
            ene.append(en)
    print(ene)
    ### matplotlib


    return 0

def main():
    parser = argparse.ArgumentParser(description='To get rid of abnormal DOS at start energy')
    parser.add_argument('inf')
    args = parser.parse_args()
    
    plot_gibbs(args.inf) 

if __name__ == '__main__':
    main()
