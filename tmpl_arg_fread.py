#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os

def modify_dos(fname, ext):
    ofname = fname+'tmp'
    fout=open(ofname, 'w')
    with open(fname, 'r') as f:
        ### analyze DOSCAR
        for i, line in enumerate(f):          # line has "\n"
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
    parser = argparse.ArgumentParser(description='To get rid of abnormal DOS at start energy')
    parser.add_argument('file', nargs='?', default='DOSCAR', help='read DOSCAR')
    parser.add_argument('-o', '--old', default='_o', help='save original DOSCAR to DOSCAR_o')
    args = parser.parse_args()
    
    modify_doscar(args.file, args.old) 

if __name__ == '__main__':
    main()
