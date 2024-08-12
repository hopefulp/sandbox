#!/home/joonho/anaconda3/bin/python
### written by David Park

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

    ###necessary lists###
    rowindexes=[]
    colindexes=[]
    whatcol=[]#selection of data
    onedarray=[]#1D-array of data
    twodarray=[]#2D-array of data
    datanames=[]

    ###use pandas dataframe
    df = pd.read_csv(fname)

    ###making list of data selection
    if y_columns == []:
        total_rows=len(df.axes[0])
        total_cols=len(df.axes[1])-1
        for i in range(1,len(df.axes[1])):
            whatcol.append(i)
    else:
        total_rows=len(df.axes[0])
        whatcol=y_columns
        total_cols=len(whatcol)

    ###making rowindexes list
    for i in range(total_rows):
        rowindexes.append(df.iloc[i][0])
    
    ###making colindexes list
    for i in (whatcol):
        for j in range(0,total_rows):
            colindexes.append(df.iloc[j][i])
    
    ###data names for legends label
    for i in (legends):
        datanames.append(i)


    ###making list of color selection###
    colorselect=[]
    if colors==[]:
        for i in range(len(whatcol)):
            colorselect.append('k')
    else:
        for i in colors:
            colorselect.append(i)

    ###making data from 1D-array to 2D-array###
    ###twodarray = 2D-array of data###
    onedarray=np.asarray(colindexes)
    twodarray=onedarray.reshape((total_cols,total_rows))
    

    ###error if number of color and data mismatches
    if len(colorselect)!=total_cols:
        print('Error: Please match the number of data and color')
        sys.exit(1)

    ###matplotlib###
    fig = plt.figure()
    ax = fig.add_subplot()
    plt.ylabel(ylabel,rotation=90, labelpad=18)
    plt.xlabel(xlabel)
    plt.title(title)


    ###dataline draw###
    if len(datanames) > 0:
        for i in range(0,total_rows):
            for j in range(0,total_cols):
                if i == 0:
                    plt.plot([i*2 ,i*2+1], [twodarray[j][i],twodarray[j][i]],color=colorselect[j],label='{0}'.format(datanames[j]))
                else:
                    plt.plot([i*2 ,i*2+1], [twodarray[j][i],twodarray[j][i]],color=colorselect[j])
    if len(datanames) == 0:
        for i in range(0,total_rows):
            for j in range(0,total_cols):
                plt.plot([i*2 ,i*2+1], [twodarray[j][i],twodarray[j][i]],color=colorselect[j])

    ###connection line draw###
    for i in range(0,total_rows-1):
        for j in range(0,total_cols):
            plt.plot([i*2+1,i*2+2],[twodarray[j][i],twodarray[j][i+1]],color=colorselect[j],linestyle='--',linewidth=0.5)


    if legends==[]:
        pass
    else:
        plt.legend()





    ax.plot()
    
    ###xticks setting(x axis data names)
    xx=list(np.arange(0.5,total_rows*2,2))
    ax.set_xticks(xx)
    ax.set_xticklabels(rowindexes, fontsize=8)
    

    plt.show()
    return 0

def main():
    parser = argparse.ArgumentParser(description='To get rid of abnormal DOS at start energy')
    parser.add_argument('inf')
    gplot = parser.add_argument_group(title='plot')
    gplot.add_argument('-yi','--y_columns',required=False, default=[], type=int, nargs='+', help='select y column index')
    gplot.add_argument('-xl', '--xlabel', default='', type=str, help='input x-label')
    gplot.add_argument('-yl', '--ylabel', required=False, default='G [eV]', help='input y-label')
    gplot.add_argument('-l', '--legends', required=False, type=str, default=[], nargs='+', help='input legends')
    gplot.add_argument('-t', '--title', required=False, type=str, default='', help='input title')
    gplot.add_argument('-c', '--color', default=['r','m','b','g'], nargs='+', help='input color')

    args = parser.parse_args()
    
    plot_gibbs(args.inf, args.y_columns, args.xlabel, args.ylabel, args.legends, args.title, args.color) 

if __name__ == '__main__':
    main()
