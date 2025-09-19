#!/home/joonho/anaconda3/bin/python
'''
    2021.02.22 pandas is imported to sort x-values with ccdd
'''
import argparse
import re
import sys
import numpy as np
from mplot2D import mplot_nvector
from libfile import load_table

def plot_pd(fin, icx, jcy, plot_dict, sep=None):

    # check inputs
    #print(f"jcy indices {jcy}")

    x, xlabel, ys, yls = load_table(fin, icx, jcy, sep)

    if not 'legends' in plot_dict.keys():
        plot_dict['legends'] = yls
    if not 'xlabel' in plot_dict:
        plot_dict['xlabel'] = xlabel
    
    mplot_nvector(x, ys, plot_dict=plot_dict, Lsave=True, v_legend=yls)

    return 0 


def main():
    parser  = argparse.ArgumentParser(description='line extraction from file and plot')
    parser.add_argument( 'fname', help='input file with tabular data')
    parser.add_argument( '-i', '--icx', default=0, type=int, help='index for x-axis, starting from 0')
    parser.add_argument( '-j', '--jcy', default=[1], type=int, nargs='+', help='set column in input file')
    parser.add_argument( '-s', '--sep', help='set column in input file')
    #parser.add_argument( '-s', '--save', action="store_true", help='save figure option')
    dfplot    = parser.add_argument_group(title = 'pandas plot')
    dfplot.add_argument( '-t', '--title', default='DFT-alpha', help='title')
    dfplot.add_argument('-pt',  '--plot_type', default='plot', choices=['plot','scatter','histo'], help="plot type") 
    dfplot.add_argument('-xs', '--xspacing', type=int, default=2, help="x-ticks")
    dfplot.add_argument('-xl', '--xlabel', default='r (A)', help='xlabel for diss distance in A')
    dfplot.add_argument('-yl', '--ylabel', default='E [eV]', help='ylabel for energy (eV)')
    dfplot.add_argument('-c', '--colors', nargs='*', default=['r','g','b','k'], help='colors')
    dfplot.add_argument('-xi', '--xlim', nargs=2, type=float, help='xrange xmin, xmax')
    dfplot.add_argument('-yi', '--ylim', nargs=2, type=float, help='yrange ymin, ymax')
    dfplot.add_argument('-xst', '--xs_type', default='numerical', choices=['evenly','numerical'], help="x-ticks spacing rule")
    dfplot.add_argument( '-dl', '--datalabel', help='how to write data labels, columns name in df for only y-values')

    args = parser.parse_args()
    mpl_xticks={}
    mpl_xticks['xspacing'] = args.xspacing
    mpl_xticks['xs_type'] = args.xs_type

    plot_dict={}
    if args.xlabel: plot_dict['xlabel'] = args.xlabel
    if args.ylabel: plot_dict['ylabel'] = args.ylabel
    if args.xlim:   plot_dict['xlim']   = args.xlim
    if args.ylim:   plot_dict['ylim']   = args.ylim
    if args.title:  plot_dict['title']  = args.title
    if args.colors: plot_dict['colors'] = args.colors
    if args.datalabel: plot_dict['legends'] = args.datalabel
    if args.plot_type:  plot_dict['type']   = args.plot_type
    if 'mpl_xticks' in locals(): plot_dict['xtick']  = mpl_xticks

    
    ### use pandas or not
    plot_pd(args.fname, args.icx, args.jcy, plot_dict, sep=args.sep)

    return 0

if __name__ == "__main__":
    main()
