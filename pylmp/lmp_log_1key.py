#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import find_file
from myplot2D import plot_line

find_kw1 = "Step"

def analyze_kw(logfile, kw):
    with open(logfile, 'r') as f:
        lines = f.readlines()
        lstep=[]
        kvalue=[]
        tag_find='OFF'
        for line in lines:
            if tag_find == 'OFF' and re.search(find_kw1,line):
                line.strip()
                lline = line.split()
                if lline[1] == find_kw1:
                    lstep.append(int(lline[2]))
                else:
                    print(f"index error in {line}")
                #print(line)
                tag_find = 'ON'
            elif tag_find == 'ON' and re.search(kw, line):
                line.strip()
                lline = line.split()
                if len(kvalue) == 0:
                    ind = lline.index(kw)
                kvalue.append(float(lline[ind+2]))
                tag_find = 'OFF'
                #print(line)
    ### check values: timestep value                
    #for x, y in zip(lstep, kvalue):
    #    print(x, y)
    ### plot 2D
    plot_line(lstep, kvalue, Title=kw, Xtitle=find_kw1, Ytitle=kw)

    return 0

def main():
    parser = argparse.ArgumentParser(description="Analyze the log.lammps")
    parser.add_argument('-l','--logfile', help="input lammps logfile ")
    parser.add_argument('-k','--kw', default="Temp", help="keyword to analyze")
    args = parser.parse_args()

    pwd = os.getcwd()
    if args.logfile == None:
        fname = find_file(pwd, 'log')
        print(f"read {fname}")
    else:
        fname = args.logfile

    analyze_kw(fname, args.kw)

if __name__ == "__main__":
    main()
