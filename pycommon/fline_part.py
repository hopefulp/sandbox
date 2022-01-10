#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os
import sys
from common import *
from f_kw import get_kw

def f_lines_extract(fname, jobtype, kw1, kw2, fext, opts, Lrun):
    tag=0
    f_name = re.split('\.', fname)
    newfile=f_name[0]+'.'+fext

    if jobtype == 'band':
        ind = opts[0]
        newfile += f'{ind}'
        if len(opts) == 2:
            fnd = opts[1]
            newfile += f'-{fnd}'
    print(f"output file is {newfile}")
        
    ### extract exact key words
    fw=open(newfile, 'w')
    ### how to read part
    with open(fname, 'r') as f:
        for line in f:
            #print(f'tag {tag}: line {line}')
            if tag == 0:
                if not kw1 in line:
                    #print(f'kw1 {kw1} not in : {line.strip()}')
                    continue
                ### if kw1 in line
                else:
                    #print(f'tag {tag}: line {line.strip()}')
                    if jobtype == 'band':
                        bandind = line.strip().split()[-1]
                        #print(f"bandline: {bandline}")
                        if bandind == ind:
                            tag = 1
                            print(line, end='')
                            fw.write(line)
                        continue
                    else:
                        tag = 1
                        continue
            #### if tag == 1: write to newfile 
            else:
                print(f'tag {tag}: line {line.strip()}')
                if jobtype == 'band':
                    if line.strip():
                        print(line, end='')
                        fw.write(line)
                    else:
                        break   # extract only one band
                else:
                    if not kw2 in line:
                        print(f"{line}, end=''")
                        fw.write(line)
                    #### if tag == 'ON' and kw is in line
                    else:
                        print(f"{line}")
                        fw.write(line)
                        break
    fw.close()
    return 0

def main():
    parser = argparse.ArgumentParser(description='file management: read file, write new file')
    parser.add_argument('fname', help='get input filename')
    parser.add_argument( '-j', '--jobtype', choices=['molden','band'],  help='job types to read a part of a file')
    parser.add_argument( '-k1', '--key1', help='using 1st key-word or "qcmolden" for molden input')
    parser.add_argument( '-k2', '--key2', help='using 2nd key-word')
    xmgroup = parser.add_argument_group(title = "XMgrace")
    xmgroup.add_argument('-i', '--index', nargs='+', default=['25'], help='extract parts using band index: start index')
    parser.add_argument( '-o', '--ext', help='output filename extension')
    parser.add_argument( '-r', '--run', action='store_true', help='run or not-False')
    args = parser.parse_args()
    
    ### define Job type
    if args.jobtype:
        jobtype = args.jobtype
    else:
        if re.search('BAND', args.fname, re.I):
            jobtype = 'band'
        else:
            jobtype = 'molden'
    print(f"jobtype: {jobtype}")            
    ### in case start & end keywords are different
    
    if args.key1:
        kw1 = args.key1
        if args.key2:
            kw2 = args.key2
        else:
            kw2 = ''
    else:
        kw1, kw2 = get_kw(jobtype)

    #if not 'kw2' in locals():
    #    kw2 = kw1
    #    print(f"kw2 == kw1 {kw1}")

    print(f"keyword 1: {kw1}")
    print(f"keyword 2: {kw2}")

    if args.ext:
        ext = args.ext
    else:
        if jobtype == 'molden':
            ext = 'molden'
        else:
            ext = 'ibd'
    f_lines_extract(args.fname, jobtype, kw1, kw2, ext, args.index, args.run)

if __name__ == "__main__":
    main()
