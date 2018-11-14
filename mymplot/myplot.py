#!/gpfs/home/joonho/anaconda3/bin/python

import argparse
import re
from _mplot2d import _mplot_2f3c

def f_sub(files,ncol,title,xt,yt,Lsave):
    nfile = len(files)
    xl  = []      # x for time
    y1 = []     # y for energy
    y2 = []

    #print("%d %d" % (nfile, ncol))
    # as for f1(x,y), f2(x,y), plot x, f1(y), f2(y)
    if nfile == 2 and ncol == 2:
        with open(files[0],"r") as f, open(files[1],"r") as g:
            for line1, line2 in zip(f, g):
                line1.strip()
                line2.strip()
                l_line1 = line1.split()
                l_line2 = line2.split()
                xl.append(float(l_line1[0]))
                y1.append(float(l_line1[1]))
                y2.append(float(l_line2[1]))
    #for x, y, z in zip(xl, y1, y2):
    #    print("%10.5f%15.5f%15.5f" % (x, y, z))

    _mplot_2f3c(xl, y1, y2, files[0], files[1],title,title,xt,yt, Lsave)

    return 0

def main():
    parser = argparse.ArgumentParser(description='2D plot w. list of files "-f f1 f2 f3 ..."')
    #parser.add_argument('job', choices=["ene"], help='type of input for gromacs')
    parser.add_argument('-f', '--files', nargs='+',help='add all the files')
    parser.add_argument('-nc', '--ncolumn', default=2, type=int, help='number of columns in each file')
    parser.add_argument('-t', '--title', default='test', help='title of figure == filename')
    parser.add_argument('-x', '--xtitle', default='t(ps)', help='X title')
    parser.add_argument('-y', '--ytitle', default='E(kJ/mol)', help='Y title')
    parser.add_argument('-s', '--save', action='store_true', help='Save figure')
    args = parser.parse_args()
   
    #print(args.files)
    f_sub(args.files,args.ncolumn,args.title,args.xtitle,args.ytitle,args.save)
    return 0

if __name__=='__main__':
	main()	


