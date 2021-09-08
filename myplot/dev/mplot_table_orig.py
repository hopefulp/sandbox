#!/home/joonho/anaconda3/bin/python
''' read file and plot '''

import argparse
import re
import sys
import numpy as np
from myplot2D import mplot_levels
from plot_job import get_jobtitle
from plot_job import Ni_6x
from my_chem import *
from common import *
import parsing 


def convert2value(st):
    try:
        val = float(st)
    except:
        if re.search('j', st):
            val = j2cal
            if re.search('-', st):
                val *= -1
        else:
            print("Error:: No transform unit for y-scale")
            sys.exit(2)
    #print(f"return scale {val}")
    return val


def get_yv_scale(yscale):
    ''' treat yscale, values as list '''
    values=[]
    for ys in yscale:
        val = convert2value(ys)
        values.append(val)
    return values

def get_title(name):
    title = fname_pre(name)
    return title.upper()

"""
    This works for one file with several y-values
    if x is included in file, modify
    if 1st line is label, modify
"""    
def draw_table(inf,fmt,icx,icy,job,title,xlabel,ylabel,line_label,Lsave,yscale,colors,Ltwinx,icy_right,ylabel_r):
    '''
    flist: file list
    '''
    ### modify get title
    if not title:
        title = inf.split('.')[0]
    
    x=[]
    y2d=[]
    ### scan file list
    with open(inf,"r") as f:
        lines=f.readlines()
        table_lines = []     # y as 2d [ [y1], [y2], ...] this needs to be transposed
        ### if 1 file, get line-labels here
        for line in lines:
            items = line.strip().split()
            table_lines.append(items)
    for i, tabline in enumerate(table_lines):
        if i==0:
            ylegend = tabline
        else:
            #print(f"use icx {icx}")
            x.append(tabline.pop(0))
            y2d.append(tabline)

    y2 = np.array(y2d).T

    print(f"size: x {len(x)}, y: {len(ylegend)}, shape of data {y2.shape}")
    ### var: x, y2, ylegend
    ### change string to value
    print(f"{y2}")
    y2value = [ [ None  if y.isalpha() else np.float(y) for y in ys ] for ys in y2 ]
    print(f"{y2value}")
    yplot = []
    y_legend = []
    for i in icy:
        yplot.append(y2value[i][:])
        y_legend.append(ylegend[i])
        
    print(f"title = {title} xlabel = {xlabel}")
    print(f"x {x}, yplot {yplot}, legend {y_legend}")
    ### x, yplot, y_legend are lists
    ### x will be used for just len(x)
    mplot_levels(x, yplot, title=title, xlabel=xlabel, ylabel=ylabel, legend=y_legend,Colors=colors)
    return 0

def main():
    parser = argparse.ArgumentParser(description='Drawing files of table')

    parser.add_argument('inf', help='read table from file')
    parser.add_argument('fmt', default='white', nargs='?', choices=['white', 'csv'], help='format of input file')
    parser.add_argument('-icx', '--icolumn_x', type=int, default=0, help='column index of X')
    g_file=parser.add_argument_group('Files', description="get input files")
    g_file.add_argument('-icy', '--icolumn_y', default=[0], nargs="*", type=int, help='column index of Y')
    g_file.add_argument('-ys', '--y_scale', default=[1], nargs="+", help='scale factor for Y [value|str|str-], use for str- for "-"')
    #g_file.add_argument('-ys', '--y_scale', nargs="*", help='scale factor for Y [value|str|str-], use for str- for "-"')
    g_twin = parser.add_argument_group('Twin-X', description='to plot using two y-axes')
    g_twin.add_argument('-tx', '--twinx', action="store_true", help='using two y-axes with twin x ticks')
    g_twin.add_argument('-icy2', '--second_iy', default=[2], nargs="+", type=int, help='designate the index of y for 2nd y-axis')
    g_twin.add_argument('-yl2', '--second_yl', default='E(ev)', help='input left y-axis title')
    parser.add_argument('-j', '--job', help='job of qcmo|ai|gromacs')
    parser.add_argument('-t', '--title', default='H2-adsorption on Pt-C54', help='title of figure would be filename')
    parser.add_argument('-xl', '--xlabel', default='Reaction Coordinate', help='X title, label in mpl')
    parser.add_argument('-yl', '--ylabel', default='E(eV)', help='Y title, label in mpl')
    parser.add_argument('-yls', '--ylabels', nargs='*', help='Y labels for legend')
    parser.add_argument('-c', '--colors', nargs='*', help='Y label for legend')
    parser.add_argument('-s', '--save', action='store_true', help='Save figure')
    args = parser.parse_args()

    ### n columns in 1 file, twinx
    draw_table(args.inf,args.fmt,args.icolumn_x,args.icolumn_y,args.job,args.title,args.xlabel,args.ylabel,args.ylabels,args.save,args.y_scale, args.colors, args.twinx, args.second_iy, args.second_yl)
    return 0

if __name__=='__main__':
	main()	


