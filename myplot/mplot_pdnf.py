#!/home/joonho/anaconda3/bin/python
'''
    2021.02.22 pandas is imported to sort x-values with ccdd
'''
import argparse
import re
import sys
import numpy as np
import _mplot2d as myplot
from _mplot2d import ax as ax
from _mplot2d import fig as fig
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
from common import *

def fplot(fin, icx, icy, ptype, title, mpl_xticks ):

    print(f"{ptype}: (x={icx}, y={icy})")
    x=[]
    y2d=[]
    #y_x=[]
    with open(fin, 'r') as f:
        #xy2d = []
        lines = f.readlines()
        ### scan data file
        for line in lines:
            cols = line.strip().split() # eles=='elements'
            #if len(ics)==2:
            if icx == 0:
                x.append(cols[icx])
            yrow=[]
            ### for any dim(y)
            for iy in icy:
                yrow.append(float(cols[iy]))
            y2d.append(yrow)
    ### x in array, y in 2d [ [y00, y01], [y10, y11], ... ] ] - get y[:,0], y[:,1], etc
    #ymin=min(y1)
    #ymax=max(y1)
    #y_min = int(ymin)
    #y_max = int(ymax) + 1
    #nbin = (y_max-1) * 10
    #my_mplot2d.draw_histogram(y1, nbin, Lsave, fig_file)
    y2d=np.array(y2d)
    print(f"shape of x, y2d : {np.array(x).shape} {y2d.shape}")
    #x=xy[:,0]
    y=y2d[:,0]
    ### mpls
    if 'xspacing' in mpl_xticks.keys() and mpl_xticks['xspacing'] != 1:
        myplot.ax_f(mpl_xticks)    
    if len(icy) == 1:
        print(f"dim of x, y = {np.array(x).shape} {y.shape}")
        if ptype == "histo":
            plt.hist2d(x, y, bins=200)
        else:
            plt.plot(x, y)
    #cb = plt.colorbar()
    #cb.set_label('counts in bin')
    plt.show()
    return 0 


def nfplot_pd(fin, icx, jcy, ptype, title, mpl_xticks, data_label):
    nfile=len(fin)
    dlabels=[]
    for i in range(nfile):
        ### PANDAS
        ### read_csv can read any type of delimiter
        data_label = f_root(fin[i])
        if 'dfm' in locals() and data_label in dfm.columns:
            data_label += '%s' % str(i)
        dlabels.append(data_label)
        df = pd.read_csv(fin[i], delim_whitespace=True, names=["X", data_label])
        #dfs.append(df)
        if i == 0:
            dfm = df
        else:
            dfm = pd.merge(dfm, df, on='X', how='outer')
        print(dfm)

    ### extract digits in 'X' in case X = \w+\d+
    dfm["N"] = dfm.X.str.replace(r"\D", "").astype(int)
    dfm = dfm.sort_values(by="N")
    print(dfm)

    #my_mplot2d.draw_histogram(y1, nbin, Lsave, fig_file)
    ### mpls
    #print(f"xticks spacing: {mpl_xticks['xspacing']}")
    ### categorical plotting
    #for col in dfm.columns[1:nfile+1]:
    for col in dlabels:
        dfm.plot(x="N", y=col, style='o-', ax=ax)
    ax.set_xticks(dfm.N)
    ax.set_ylabel('Frmse (eV/A)')
    ax.set_xticklabels(dfm.X)
    ax.set_title(title)
    if 'xspacing' in mpl_xticks.keys() and mpl_xticks['xspacing'] != 1:
        print(f"control xticks spacing: {mpl_xticks['xspacing']}")
        ### only xlabel.set_visible() works for Numerically-spaced
        #ax.xaxis.set_major_locator(ticker.MaxNLocator(mpl_xticks['xspacing']))
        #ax.xaxis.set_major_locator(ticker.MultipleLocator(mpl_xticks['xspacing']))
        #plt.locator_params(nbins=mpl_xticks['xspacing'])
        myplot.xtick_invisible(ax, mpl_xticks['xspacing'])

    #cb = plt.colorbar()
    #cb.set_label('counts in bin')
    plt.show()
    return 0 


def main():
    parser  = argparse.ArgumentParser(description='line extraction from file and plot')
    parser.add_argument( 'fname', nargs='+', help='input file')
    parser.add_argument( '-i', '--icx', default=0, help='index for x-axis, starting from 0')
    parser.add_argument( '-j', '--jcy', default=[1], type=int, nargs='+', help='set column in input file')
    #parser.add_argument( '-s', '--save', action="store_true", help='save figure option')
    parser.add_argument( '-pd', '--pandas', action='store_true', help="use pandas to extract number and to use ordering")
    dfplot    = parser.add_argument_group(title = 'pandas plot')
    dfplot.add_argument( '-t', '--title', default='Train-Test', help='title')
    dfplot.add_argument('-pt',  '--plot_type', default='plot', choices=['plot','scatter','histo'], help="plot type") 
    dfplot.add_argument('-xs', '--xspacing', type=int, default=2, help="x-ticks")
    dfplot.add_argument('-xst', '--xs_type', default='numerical', choices=['evenly','numerical'], help="x-ticks spacing rule")
    dfplot.add_argument( '-dl', '--datalabel', help='how to write data labels, columns name in df')

    args = parser.parse_args()
    mpl_xticks={}
    mpl_xticks['xspacing'] = args.xspacing
    mpl_xticks['xs_type'] = args.xs_type
    
    ### use pandas or not
    if args.pandas:
        nfplot_pd(args.fname, args.icx, args.jcy, args.plot_type, args.title, mpl_xticks, args.datalabel)
    else:
        fplot(args.fname, args.icx, args.icy, args.plot_type, args.title, mpl_xticks)
    return 0

if __name__ == "__main__":
    main()
