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

def fplot(fin, icx, icy, fout, mini, maxi, i_line, title, ptype, mpl_xticks ):

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


def fplot_pd(fin, icx, icy, fout, mini, maxi, i_line, title, ptype, mpl_xticks ):
    nfile=len(fin)
    f = fin[0]
    data_label=f_root(f)
    ### PANDAS
    ### read_csv can read any type of delimiter
    df = pd.read_csv(f, delim_whitespace=True, names=["X", "Y"])
    ### extract digits in 'X' in case X = \w+\d+
    df["N"] = df.X.str.replace(r"\D", "").astype(int)
    df = df.sort_values(by="N")
    ### 
    if nfile == 2:
        f2 = fin[1]
        data_label2 = f_root(f2)
        df2 = pd.read_csv(f2, delim_whitespace=True, names=["X", "Y2"])
        df = pd.merge(df, df2, on='X')
    print(df)
    #my_mplot2d.draw_histogram(y1, nbin, Lsave, fig_file)
    ### mpls
    #print(f"xticks spacing: {mpl_xticks['xspacing']}")
    if len(icy) == 1:
        ### categorical plotting
        if mpl_xticks['xs_type'] == 'evenly':
            df.plot(x="X", y="Y", ax=ax)
            ax.set_title("Evenly spaced")
        else:
            if nfile ==1:
                df.plot(x="N", y="Y", style='o-', ax=ax, label=data_label)
            else:
                ### change columns name and plot
                df2 = df.rename(columns={'Y':data_label, 'Y2':data_label2})
                df2.plot(x="N", y=[data_label,data_label2], style=['o-','o-'], ax=ax) # two labels are not working
            ax.set_xticks(df.N)
            ax.set_ylabel('Frmse (eV/A)')
            ax.set_xticklabels(df.X)
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
    parser.add_argument( '-ix', '--icx', default=0, help='index for x-axis')
    parser.add_argument( '-iy', '--icy', default=[1], type=int, nargs='+', help='set column in input file')
    parser.add_argument( '-m', '--mini', help='minimum value of random variable')
    parser.add_argument( '-n', '--maxi', help='maximum value of random variable')
    parser.add_argument( '-i', '--init', default=0, type=int, help='start line number')
    parser.add_argument( '-f', '--final', type=int, help='end line number')
    #parser.add_argument( '-s', '--save', action="store_true", help='save figure option')
    parser.add_argument( '-pd', '--pandas', action='store_false', help="use pandas to extract number and to use ordering")
    parser.add_argument( '-t', '--title', default='Train-Test', help='title')
    mpls    = parser.add_argument_group(title = 'matplotlib args')
    mpls.add_argument('-pt',  '--plot_type', default='plot', choices=['plot','scatter','histo'], help="plot type") 
    mpls.add_argument('-xs', '--xspacing', type=int, default=2, help="x-ticks")
    mpls.add_argument('-xst', '--xs_type', default='numerical', choices=['evenly','numerical'], help="x-ticks spacing rule")

    args = parser.parse_args()
    mpl_xticks={}
    mpl_xticks['xspacing'] = args.xspacing
    mpl_xticks['xs_type'] = args.xs_type
    
    ### use pandas or not
    if args.pandas:
        fplot_pd(args.fname, args.icx, args.icy, args.mini, args.maxi, args.init, args.final, args.title, args.plot_type, mpl_xticks)
    else:
        fplot(args.fname, args.icx, args.icy, args.mini, args.maxi, args.init, args.final, args.title, args.plot_type, mpl_xticks)
    return 0

if __name__ == "__main__":
    main()
