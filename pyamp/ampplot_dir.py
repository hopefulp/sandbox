#!/home/joonho/anaconda3/bin/python

import argparse
import myplot2D                ### FOR MLET: QXcbConnection: Could not connect to display
from common import f_parsing
import sys
import os
import py_method
import glob
import re
import subprocess
import amp_ini
import pickle

def get_x(xin, xv):
    if re.search('d', xin):    # for dir
        dirs = os.listdir()
    elif re.search('p', xin):
        dirs = glob.glob(f"{xv}*")
    dirs.sort(key=py_method.natural_keys)
    x = dirs
    return x

def get_y(yin, xv, sub_dir=None):
    y = []
    if re.search('amp-log', yin):
        for x in xv:
            if sub_dir:
                com = "grep 'optimization ' %s/%s/amp-log.txt -B 1 | awk '/T/ {print $8, $10}'" % (x, sub_dir)
            else:
                com = "grep 'optimization ' %s/amp-log.txt -B 1 | awk '/T/ {print $8, $10}'" % x
            results = py_method.get_cli(com)
            yv = float(results.split()[0])
            y.append(yv)
    return y

def get_y2d(yfile, xv, y_dirs):
    y2d = []
    for di in y_dirs:
        y=[]
        for x in xv:
            if yfile == 'amp-log.txt':
                subdir = f'{x}/{di}'
                if not os.path.isdir(subdir):
                    print(f"dirname {di} does not exist")
                    sys.exit(10)
                com = "grep 'optimization ' %s/%s/amp-log.txt -B 1 | awk '/T/ {print $8, $10}'" % (x, di)
                results = py_method.get_cli(com)
                yv = float(results.split()[0])
                y.append(yv)
            elif yfile == 'test_fstat_acc.pkl':
                fname = '%s/%s/test_fstat_acc.pkl' % (x, di)
                with open(fname, 'rb') as f:
                    p = pickle.load(f)
                    frmse = p.rmse
                    y.append(frmse)
        y2d.append(y)
        
    return y2d

def get_legend(legends):
    legend=[]
    for leg in legends:
        if leg == '.':
            legend.append('Ndata100')
        else:
            legend.append(leg)
    return legend

def mpyplot(xin, xv, yv, y_dirs, xlabel, ylabel, title):
    ### obtain X's
    x_orig = get_x(xin, xv)
    nx = len(x_orig)

    ### obtain Y's
    #y = get_y(yv, x_orig, y_dirs)
    y2d = get_y2d(yv, x_orig, y_dirs)
    leg_tag = y_dirs
    legend = get_legend(leg_tag)
    print(f"{x_orig}")
    print(f"{y2d}")

    ### change x in integer
    x=[]
    for xch in x_orig:
        x.append(int(re.sub("\D", "", xch)))
    myplot2D.mplot_nvector(x, y2d, title = title, xlabel=xlabel, ylabel=ylabel, legend=legend)
    return 0

def main():
    parser = argparse.ArgumentParser(description='plot amp test using my2dplot')
    xgroup = parser.add_mutually_exclusive_group(required=True)
    xgroup.add_argument('-x', '--xvalues',  help='scan present directory ')
    xgroup.add_argument('-p', '--prefix',  help='scan present directory using prefix')
    parser.add_argument('-yd', '--ydirs', nargs='+', default=['.'], help='read data in subdirectories as much as num of yvalues')
    #parser.add_argument('-yf', '--yfiles', default='amp-log.txt', choices=['amp-log.txt','test_fstat_acc.pkl'], help='amp pot out file ')
    parser.add_argument('-y', '--yvalue', default='tr', choices=['tr','te'], help='amp pot out file ')
    parser.add_argument('-xl', '--xlabel', help='xlabel for matplot')
    parser.add_argument('-yl', '--ylabel', help='ylabel for matplot')
    parser.add_argument('-tl', '--title',  help='title for matplot')
    args = parser.parse_args()

    if args.xvalues:
        xin = 'dir'
        xvalue = args.xvalues
    elif args.prefix:
        xin = 'p'
        xvalues = args.prefix

    if not args.xlabel:
        args.xlabel = "Nparam"
    if not args.ylabel:
        args.ylabel = r"F$_{rmse}$ (eV/$\AA$)"
    if not args.title:
        args.title = "F$_{rmse}$ vs Nparam (Gs-pow)"
    if args.yvalue == 'tr':
        yvalue = 'amp-log.txt'
    elif args.yvalue == 'te':
        yvalue = 'test_fstat_acc.pkl'
    else:
        yvalue = args.yvalue

    mpyplot(xin, xvalues, yvalue,args.ydirs,args.xlabel,args.ylabel,args.title)
    return 0

if __name__ == '__main__':
    main()

