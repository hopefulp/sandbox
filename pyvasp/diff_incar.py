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
        if fname != 'INCAR' and os.path.isdir(fname):
            fname = fname + '/INCAR'
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
    return dic_list

def print_key2(dic, params):
    for p in params: 
        if p.upper() in dic.keys():
            print(f"{p.upper()} {dic[p.upper()]}") 
        else: 
            print(f"{p.upper()} does not exist")
    return 0

def print_key(dic, params):
    for p in params:
        tag_search = False
        for key in dic.keys():
            if re.search(p.upper(), key):
                print(f"{key} {dic[key]}")
                tag_search = True
            if tag_search:
                break
        if tag_search:
            break
        else:
            print(f"{p.upper()} does not exist")
    return 0

def compare_incar(files, key, Ldiff=False):
    tup = get_incar_dic(files)
    if len(tup) == 1:
        if key:
            print_key(tup[0], key)
        else:
            print_dict(tup[0])
        return 0
    elif len(tup) == 2:
        d1 = tup[0]
        d2 = tup[1]
        if key:
            for t in tup:
                print_key(t)
            return 0

    d1keys = set(d1.keys())
    d2keys = set(d2.keys())
    shared = d1keys.intersection(d2keys)
    diff_1 = d1keys.difference(d2keys)
    diff_2 = d2keys.difference(d1keys)
    union  = d1keys.union(d2keys)
    s = ' '
    for key in union:
        if key in shared:
            if d1[key] == d2[key]:
                diff="-"
            else:
                diff="diff"
            print(f'{s:4}{key:<8} {d1[key]:>10} {d2[key]:>10} {diff:>10}')
        elif key in diff_1:
            diff="L"    # "1st"  
            print(f'{key:<12} {d1[key]:>10} {s:10} {s:6}{diff:<4}')
        else:
            diff="R"    # 'last'
            print(f'{key:>12} {s:10} {d2[key]:>10} {diff:>10}')

    return 0
    

def main():

    parser = argparse.ArgumentParser(description="Compare two INCAR")
    parser.add_argument('files', nargs='*',  help="one or two INCAR file ")
    parser.add_argument('-p', '--params', nargs='*',  help="check the key and values ")
    args = parser.parse_args()

    if not args.files:
        files=['INCAR']
    else:
        files=args.files

    compare_incar(files, args.params)

if __name__ == "__main__":
    main()
