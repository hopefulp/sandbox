#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import csv
import sys

def plot_gibbs(fname):
    x=[]
    ene=[]
    with open(fname, 'r') as f:
        lines=f.readlines()
        for line in lines:
            en = float(line.strip())
            ene.append(en)
#        reader = csv.DictReader(f)
#        for row in reader:
#            x.append(row['x'])
#            ene.append(float(row['y']))

    ### matplotlib
    fig=plt.figure()
    ax = plt.axes()
    plt.ylabel('E (eV)',rotation=0, labelpad=18)
    varnum=len(ene)

    for i in range(varnum):
        plt.plot([i*2 -1 ,i*2], [ene[i],ene[i]],color='k')
    for i in range(varnum -1):
        plt.plot([i*2,i*2+1],[ene[i],ene[i+1]],color='k',linestyle='--',linewidth=0.5)

    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

    plt.show()
    return 0

def main():
    parser = argparse.ArgumentParser(description='To get rid of abnormal DOS at start energy')
    parser.add_argument('inf')
    args = parser.parse_args()
    
    plot_gibbs(args.inf) 

if __name__ == '__main__':
    main()
