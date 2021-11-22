#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_files

def plot_band(inf, ymax, ymin):
    '''
    f has many lines
    '''
    kw = 'Band-Index'
    kp=[]
    eband=[]
    iband=[]
    with open(inf, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if kw in line:
            ele = line.split()
            iband.append(ele[2])

        

    if job == None:
        print("List this directory = ")
    return 0 
    

def main():

    parser = argparse.ArgumentParser(description="read many tables in a file")
    parser.add_argument('-d','--dat', default="BAND.dat", help="plot BAND.dat from vaspkit ")
    parser.add_argument('-ym','--ymax', default=2.0, type=float,  help="ylim max")
    parser.add_argument('-yn','--ymin', default=-2.0, type=float,  help="ylim min")
    args = parser.parse_args()

    plot_band(args.dat, args.ymax, args.ymin)

if __name__ == "__main__":
    main()
