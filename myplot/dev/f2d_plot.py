#!/home/joonho/anaconda3/bin/python

import argparse
import re
import sys
import numpy as np
import _mplot2d as myplot
import matplotlib.pyplot as plt

#def num(arg):
#    if arg

def fplot(fin, icx, icy, fout, mini, maxi, i_line, Lsave, ptype, mpl_xticks ):
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


def main():
    parser  = argparse.ArgumentParser(description='line extraction from file and plot')
    parser.add_argument( 'fname', help='input file')
    mpls    = parser.add_argument_group(title = 'matplotlib args')
    mpls.add_argument('-t',  '--type', default='plot', choices=['plot','scatter','histo'], help="plot type") 
    mpls.add_argument('-xs', '--xspacing', type=int, default=2, help="x-ticks")
    parser.add_argument( '-ix', '--icx', default=0, help='index for x-axis')
    parser.add_argument( '-iy', '--icy', default=[1], type=int, nargs='+', help='set column in input file')
    parser.add_argument( '-m', '--mini', help='minimum value of random variable')
    parser.add_argument( '-n', '--maxi', help='maximum value of random variable')
    parser.add_argument( '-i', '--init', default=0, type=int, help='start line number')
    parser.add_argument( '-f', '--final', type=int, help='end line number')
    parser.add_argument( '-s', '--save', action="store_true", help='save figure option')
    #parser.add_argument( '-t', '--figure', default='gmx.png', help='png file name')

    args = parser.parse_args()
    mpl_xticks={}
    mpl_xticks['xspacing'] = args.xspacing

    fplot(args.fname, args.icx, args.icy, args.mini, args.maxi, args.init, args.final, args.save, args.type, mpl_xticks)
    return 0

if __name__ == "__main__":
    main()
