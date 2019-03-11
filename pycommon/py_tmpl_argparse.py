#!/usr/bin/python

import argparse
import re
import os

### scratch keywords for AMP from IM/MM input file of the form of .fin
kw1 = "energy"
kw2 = "cartesian"

def extract_data(f_input, fit_e):
    #fout=open("t.csh", 'w')
    with open(f_input, 'r') as f:
        for line in f:          # line has "\n"
            #print line,         # print writes its own "\n"
            if re.search(ke1, line):
                strn='#PBS -l walltime='+time+':00:00'+"\n"
                fout.write(strn)
            else:
                fout.write(line)    # write does not write "\n"
    #fout.close()
    #str1='qsub t.csh'
    #print str1
    #x=os.system(str1)
    return x

def main():
    parser = argparse.ArgumentParser(description='get coordinate [energy, force, hessian] from IM input file (.fin)')
    parser.add_argument('f_immm_input', help='IM/MM input file')
    parser.add_argument('-t', '--fit_type', default='e', help='data type: energy alone, energy & force, energy, force & hessian')
    args = parser.parse_args()
    extract_data(args.f_immm_input, args.fit_type) 

if __name__ == '__main__':
    main()
