#!/home/jackjack5/epd/bin/python

import argparse
import re
import os
import sys

def cmd(fqcout, f_kind, keyword, n_number, rmline):
    fname = re.split("\.", fqcout)
    if not len(fname) == 2:
        print "filename type error: it should have one period"
        sys.exit(1)
    if f_kind == 'irc':
        Turn = 2
    foutname = fname[0] + "-iq.out"
    fout=open(foutname, 'w')
    i = 0
    tag="OFF"
    ### 1 for remove last problematic sentence
    ### 0 for remove problematic sentence after finding keyword
    method=1

    if method == 1:
        print("at the moment, remove last \'%s problematic\' line" % rmline)
        with open(fqcout, 'r') as f:
            for line in f:      # line has "\n"
                if re.search(rmline, line):
                    i += 1
        ntime = i
        i = 0
        with open(fqcout, 'r') as f:
            for line in f:
                if re.search(rmline, line):
                    i += 1
                if i == ntime:
                    i = 1000
                    pass
                else:
                    fout.write(line)

        
    ### method to find convergence keyword
    elif method == 0:
        with open(fqcout, 'r') as f:
            for line in f:          # line has "\n"
                if tag == "OFF":
                    if re.search(keyword, line):
                        if i==0:
                            i += 1
                        elif i == 1:
                            tag="ON"
                    fout.write(line)
                else:
                    if re.search(rmline, line):
                        pass
                    else:
                        fout.write(line)
    fout.close()
    return 

def main():
    parser = argparse.ArgumentParser(description='remove last removable line for irc display in iqmol')
    parser.add_argument('outfile', type=str, help='qchem outfile')
    parser.add_argument('-f', '--file', type=str, default='irc', help='qchem irc file modification')
    parser.add_argument('-k', '--keyword', type=str, default='IRC -- convergence', help='word to be removed')
    parser.add_argument('-rm', '--remove', type=str, default='Standard Nuclear Orientation', help='word to be removed')
    parser.add_argument('-n', '--number', type=int, default=2, help='in case the keyword appears once, use 1')
    args = parser.parse_args()
    cmd(args.outfile, args.file, args.keyword, args.number, args.remove) 
    return 

if __name__ == '__main__':
    main()
