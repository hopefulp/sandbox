#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os
import sys

key_match = {'GMAX': 'local part'}
nlines = {'local part': 1}

def fextract(fname, kw):
    if os.path.isdir(fname):
        fname += '/POTCAR'
    kwlines=[]
    kvs=[]
    kw = kw.upper()
    if kw in key_match.keys():
        kw_orig = kw
        kw = key_match[kw]
    ### READ POTCAR
    with open(fname, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):          # line has "\n"
        #print line,         # print writes its own "\n"
        if re.search(kw, line):
            if kw in nlines.keys():
                kwlines.append(lines[i+nlines[kw]].strip())
            else:
                kwlines.append(line.strip())
            print(line.strip())
    if not kwlines:
        print(f"there is no kw {kw} in {fname}")
        return 1
    else:
        for kwline in kwlines:
            lstr = re.split(r'[\s;]+', kwline)      #r'[;\s]+')
            #print(lstr)
            ### if key = value line
            if lstr[0] == kw:
                kvs.append(lstr[2])
            ### if only value line
            else:
                kvs.append(lstr[0])
    if kvs:
        print(f"{kvs}")
        if kw == 'ENMAX':
            enmax = list(map(float, kvs))
            print(f"{enmax}")
            default = round((max(enmax) * 1.3)/50)*50
            print(f"{kw}: {max(enmax)} roundoff with 30% for cell relaxation: {default}")
        else:
            if len(kvs) == 1:
                print(f"{kw_orig}: {kvs[0]}")
            else:
                print(f"{kw_orig}: {kvs}")
    return 0

def main():
    parser = argparse.ArgumentParser(description='READ POTCAR')
    parser.add_argument('file', nargs='?', default='POTCAR', help='read POTCAR')
    parser.add_argument('-k', '--kw', default='enmax', help='extract value using key such as enmax,zval, gmax')
    parser.add_argument('-u', '--usage', action='store_true', help='print usage')
    args = parser.parse_args()

    if args.usage:
        print(f"potcar_kw.py POTCAR_H -k gmax 'only for local part'\
            \n\tdata_treat.py -i POTCAR2.H/0/POTCAR POTCAR2.H/0.1/POTCAR -l 57 25\
            \npotcar_kw.py POTCAR_H -k enmax\
            ")
        sys.exit(0)
    
    fextract(args.file, args.kw) 

if __name__ == '__main__':
    main()
