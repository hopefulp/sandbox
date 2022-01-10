#!/usr/bin/python

import argparse
import re
import os
import sys
from common import *

def f_manage(fname, kw1, kw2, fext, Lrun):
    tag='OFF'
    f_name = re.split('\.', fname)
    new_file=f_name[0]+'.'+fext
    print "output file is %s" % new_file
    fw=open(new_file, 'w')

    ### in case start & end keywords are different
    if kw1 == "qcmolden":
        import qc_out_kw
        kw1 = qc_out_kw.qcout_molden1
        kw2 = qc_out_kw.qcout_molden1b

    if not kw2:
        print "kw2 == kw1 %s" % kw1
        kw2 = kw1

    print "keyword 1: %s" % kw1,
    print "keyword 2: %s" % kw2
    with open(fname, 'r') as f:
        for line in f:
            if tag == 'OFF':
                if not kw1 in line:
                    #print line,
                    continue
                else:
                    tag = 'ON'
                    print "tag = ", tag
                    continue
            #### if tag == 'ON'                    
            elif not kw2 in line:
                print line,
                fw.write(line)
            #### if tag == 'ON' and kw is in line
            else:
                print line
                fw.write(line)
                break
    fw.close()
    return 0

def main():
    parser = argparse.ArgumentParser(description='file management: read file, write new file')
    parser.add_argument('fname', help='get input filename')
    #parser.add_argument( '-j','--job', choices=['cut'],  help='cut for file cut')
    parser.add_argument( 'keyword', default='qcmolden', choices=['qcmolden'], help='using 1st key-word or "qcmolden" for molden input')
    parser.add_argument( '-k2', '--key2', help='using 2nd key-word')
    parser.add_argument( '-o', '--ext', default='molden', help='output filename extension')
    parser.add_argument( '-r', '--run', action='store_true', help='run or not-False')
    args = parser.parse_args()
    #print args.key2

    f_manage(args.fname, args.keyword, args.key2, args.ext, args.run)

if __name__ == "__main__":
    main()
