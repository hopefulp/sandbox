#!/home/joonho/anaconda3/bin/python

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def test_xtick(fin):

    li=[]
    with open(fin, 'r') as f:
        lines=f.readlines()
        for line in lines:
            l =  line.strip().split()
            li.append([l[0], l[1]])

    ar = np.array(li)

    fig, ax = plt.subplots(1,1)
    x=ar[:,0]
    y=ar[:,1]
    y = y.astype(np.float)

    ### using MultipleLocator, tick is shifted
    tick_spacing=2
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    ax.plot(x, y, 'o-')
    plt.show()
    return 0

def main():
    parser  = argparse.ArgumentParser(description='line extraction from file and plot')
    parser.add_argument( 'fname', help='input file')

    args = parser.parse_args()

    test_xtick(args.fname)
    return 0

if __name__ == '__main__':
    main()

