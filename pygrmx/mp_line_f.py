#!/gpfs/home/joonho/anaconda3/bin/python

import argparse
from mplot2d import *
def main():
    parser = argparse.ArgumentParser(description='mplot 2d column')
    parser.add_argument('fname', help='input file')
    parser.add_argument('-d', '--decimal', type=int, default=1, help='number of decimal places')
    parser.add_argument('-s', '--save', action='store_true', help='save figure as png file')
    parser.add_argument('-f', '--figname', default='line.png', help='input figure filename')
    parser.add_argument('-t', '--title', default='Figure', help='figure title')
    parser.add_argument('-x', '--xtitle', default='X', help='x-axis title')
    parser.add_argument('-y', '--ytitle', default='Y', help='y-axis title')

    args = parser.parse_args()

    f_draw(args.fname, args.decimal, args.save, args.figname, args.title, args.xtitle, args.ytitle)

    return 0

if __name__ == "__main__":
    main()
