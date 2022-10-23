#!/home/joonho/anaconda3/bin/python
### written by David Park Sep. 2022

import pandas as pd
import argparse
import re
#import os
#import matplotlib as mpl
import matplotlib.pyplot as plt
#import csv
import sys
import numpy as np
'''
    This reads only csv format
'''

def plot_gibbs(fname, y_columns, xlabel, ylabel, legends, title, colors):
    ### check input file whether it is csv format
    l_fname = re.split('\.', fname)
    if l_fname[1] != 'csv':
        print("Please input csv format")
        sys.exit(1)
    ### use pandas dataframe
    df = pd.read_csv(fname)
    total_rows=len(df.axes[0])#number of row
    total_cols=len(df.axes[1])#number of column

    rowindexes=[]
    colindexes=[]

    for i in range(total_rows):
        rowindexes.append(df.loc[i][0])

#    print(rowindexes) <--- Check Complete, Successful!

    for i in range (0,total_rows):
        for j in range (1,total_cols):
            colindexes.append(df.loc[i][j])

    print(rowindexes,colindexes,total_cols)

    '''
    xs=[]
    ene=[]
    with open(fname, 'r') as f:

#       lines=f.readlines()
#       for line in lines:
#           en = float(line.strip())
#           ene.append(en)
        reader = csv.DictReader(f)
        for row in reader:
            xs.append()
            ene.append()
    '''


    ### matplotlib 
    fig = plt.figure()
    ax = fig.add_subplot()

    plt.ylabel('E (eV)',rotation=0, labelpad=18)
    varnum=total_rows






#data line
    for i in range(varnum):
        for j in range((total_cols-1)*i,(total_cols-1)*i+(total_cols-1)):
            plt.plot([i*2 ,i*2+1], [colindexes[j],colindexes[j]],color='k')
#connection line
    for i in range(varnum-1):#0~9
        for j in range(total_cols-1):#4
            plt.plot([i*2+1,i*2+2],[colindexes[(total_cols-1)*i+j],colindexes[((total_cols-1)*(i+1))+j]],color='k',linestyle='--',linewidth=0.5)
#erase all xticks
#    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
#adding xticks
    ax.plot()

    xx=list(np.arange(0.5,varnum*2,2))
    ax.set_xticks(xx)
    ax.set_xticklabels(rowindexes, fontsize=8)
        
        




    plt.show()
    return 0

def main():
    parser = argparse.ArgumentParser(description='To get rid of abnormal DOS at start energy')
    parser.add_argument('inf')
    parser.add_argument('-yi','--y_columns', default=0, type=int, narg='+', help='select y column index')
    parser.add_argument('-xl', '--xlable', help='input x-label')
    parser.add_argument('-yl', '--ylable', default='G [eV]', help='input x-label')
    parser.add_argument('-l', '--legend', help='input x-label')
    parser.add_argument('-t', '--title', help='input x-label')
    parser.add_argument('-c', '--color', help='input x-label')

    args = parser.parse_args()
    
    plot_gibbs(args.inf, args.y_columns, args.xlabel, args.ylabel, args.legend, args.title, args.color) 

if __name__ == '__main__':
    main()
