#!/gpfs/home/joonho/anaconda3/bin/python
#/usr/bin/python

import argparse
import re
import numpy as np
import my_mplot2d

def histogram(fin, fout, mini, maxi, i_line, f_line, Lsave, fig_file):
    #x1=[]
    y1=[]
    y_x=[]
    with open(fin, 'r') as f:
        i=0
        for line in f:
            i+=1
            if i < i_line:
                continue
            if re.search("[a-zA-Z]", line):
                continue
            xy=line.split()
            if not xy[0]:
                del xy[0] 
            #x1.append(float(xy[0]))
            y1.append(float(xy[1]))

    ymin=min(y1)
    ymax=max(y1)
    #print(ymin, ymax)

    y_min = int(ymin)
    y_max = int(ymax) + 1

    nbin = (y_max-1) * 10

    mplot2d.draw_histogram(y1, nbin, Lsave, fig_file)

    return 0 


def main():
    parser = argparse.ArgumentParser(description='line extraction from file and plot')
    parser.add_argument( 'fname', help='input file')
    parser.add_argument( '-o', '--outf', default='histo.dat', help='frequency in the bin')
    parser.add_argument( '-m', '--mini', help='minimum value of random variable')
    parser.add_argument( '-n', '--maxi', help='maximum value of random variable')
    parser.add_argument( '-i', '--init', default=0, type=int, help='start line number')
    parser.add_argument( '-f', '--final', type=int, help='end line number')
    parser.add_argument( '-s', '--save', action="store_true", help='save figure option')
    parser.add_argument( '-t', '--figure', default='gmx.png', help='png file name')

    args = parser.parse_args()

    histogram(args.fname, args.outf, args.mini, args.maxi, args.init, args.final, args.save, args.figure)
    return 0

if __name__ == "__main__":
    main()
