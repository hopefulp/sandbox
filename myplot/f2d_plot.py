#!/home/joonho/anaconda3/bin/python

import argparse
import re
import sys
import numpy as np
import my_mplot2d
import matplotlib.pyplot as plt

def histogram(fin, ics, fout, mini, maxi, i_line, f_line, Lsave, fig_file):
    #x1=[]
    #y1=[]
    #y_x=[]
    with open(fin, 'r') as f:
        #i=0
        xy2d = []
        lines = f.readlines()
        for line in lines:
            col = line.split()
            if len(ics)==2:
                row = [float(col[ics[0]-1]),float(col[ics[1]-1])]
                xy2d.append(row)
            else:
                print("this draws 2 columns of a file")
                sys.exit(0)
    #ymin=min(y1)
    #ymax=max(y1)
    #y_min = int(ymin)
    #y_max = int(ymax) + 1
    #nbin = (y_max-1) * 10
    #my_mplot2d.draw_histogram(y1, nbin, Lsave, fig_file)
    xy=np.array(xy2d)
    print(f"shape of xy: {xy.shape}")
    x=xy[:,0]
    y=xy[:,1]
    print(f"dim of x, y = {x.shape} {y.shape}")
    plt.hist2d(x, y, bins=200)
    cb = plt.colorbar()
    cb.set_label('counts in bin')
    plt.show()
    return 0 


def main():
    parser = argparse.ArgumentParser(description='line extraction from file and plot')
    parser.add_argument( 'fname', help='input file')
    parser.add_argument( '-c', '--column', type=int, nargs='+', help='set column in input file')
    parser.add_argument( '-o', '--outf', default='histo.dat', help='frequency in the bin')
    parser.add_argument( '-m', '--mini', help='minimum value of random variable')
    parser.add_argument( '-n', '--maxi', help='maximum value of random variable')
    parser.add_argument( '-i', '--init', default=0, type=int, help='start line number')
    parser.add_argument( '-f', '--final', type=int, help='end line number')
    parser.add_argument( '-s', '--save', action="store_true", help='save figure option')
    parser.add_argument( '-t', '--figure', default='gmx.png', help='png file name')

    args = parser.parse_args()

    histogram(args.fname, args.column, args.outf, args.mini, args.maxi, args.init, args.final, args.save, args.figure)
    return 0

if __name__ == "__main__":
    main()
