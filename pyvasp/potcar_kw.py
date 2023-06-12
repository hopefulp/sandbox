#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os

def fextract(fname, kw):
    if os.path.isdir(fname):
        fname += '/POTCAR'
    kwlines=[]
    kvs=[]
    kw = kw.upper()
    ### READ POTCAR
    with open(fname, 'r') as f:
        for i, line in enumerate(f):          # line has "\n"
            #print line,         # print writes its own "\n"
            if re.search(kw, line):
                kwlines.append(line)
                print(line.strip())
    if not kwlines:
        print(f"there is no kw {kw} in {fname}")
        return 1
    else:
        for line in kwlines:
            lstr = re.split(r'[\s;]+', line.strip()) #r'[;\s]+')
            #print(lstr)
            if lstr[0] == kw:
                kvs.append(lstr[2])
    if kvs:
        print(f"{kvs}")
        enmax = list(map(float, kvs))
        print(f"{enmax}")
        default = round((max(enmax) * 1.3)/50)*50
        print(f"{kw}: {max(enmax)} roundoff with 30% for cell relaxation: {default}")
    return 0

def main():
    parser = argparse.ArgumentParser(description='READ POTCAR')
    parser.add_argument('file', nargs='?', default='POTCAR', help='read POTCAR')
    parser.add_argument('-k', '--kw', default='ENMAX', help='extract keyword & value')
    args = parser.parse_args()
    
    fextract(args.file, args.kw) 

if __name__ == '__main__':
    main()
