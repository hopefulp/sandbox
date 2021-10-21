#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_files
from my_print import print_dict

reg='[=\s]+'

def get_kv(strline):
    lst = re.split(reg, strline)
    return lst[0], lst[1]

def get_incar_dic(files):
    dic_list = []
    for fname in files:
        dic={}
        if os.path.isfile(fname):
            with open(fname, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    ### remove white space on both edge
                    strline=line.strip()
                    ### select a line which starts with word
                    if re.match('\w',strline):
                        if re.search('=', strline):
                            if re.search(';', strline):
                                kvslist = re.split(';', strline)
                                for kvs in kvslist:
                                    strl = kvs.strip()
                                    #print(f"kvstring {strl}")
                                    k, v = get_kv(strl)
                                    dic[k] = v
                            ### delimiter: = and \s(white space) altogether
                            else:
                                k, v = get_kv(strline)
                                #print(f'{lst[0]:^10} = {lst[1]:>10}')
                                dic[k] = v
        dic_list.append(dic)                            
    return (*dic_list, )

def compare_incar(files, Ldiff=False):
    tup = get_incar_dic(files)
    if len(tup) == 1:
        print_dict(tup[0])
        return 0
    elif len(tup) == 2:
        d1 = tup[0]
        d2 = tup[1]
    d1keys = set(d1.keys())
    d2keys = set(d2.keys())
    shared = d1keys.intersection(d2keys)
    diff_1 = d1keys.difference(d2keys)
    diff_2 = d2keys.difference(d1keys)
    union  = d1keys.union(d2keys)
    s = ' '*10
    for key in union:
        if key in shared:
            if d1[key] == d2[key]:
                diff="-"
            else:
                diff="diff"
            print(f'{"    "}{key:<8} {d1[key]:>10} {d2[key]:>10} {diff:>10}')
        elif key in diff_1:
            diff=""    # "1st"  
            s=' '*10
            print(f'{key:<12} {d1[key]:>10} {s:10} {diff:^10}')
        else:
            diff=""    # 'last'
            print(f'{key:>12} {s:10} {d2[key]:>10} {diff:>10}')

    return 0
    

def main():

    parser = argparse.ArgumentParser(description="Compare two INCAR")
    parser.add_argument('files', nargs='*',  help="one or two INCAR file ")
    args = parser.parse_args()

    if not args.files:
        files=['INCAR']
    else:
        files=args.files

    compare_incar(files)

if __name__ == "__main__":
    main()
