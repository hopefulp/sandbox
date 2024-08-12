#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_files
from my_print import print_dict
from mod_incar import ordered_incar
reg='[=\s]+'

def get_kv(strline):
    lst = re.split(reg, strline)
    return lst[0], lst[1]

def get_incar_dict(f):
    dic={}
    with open(f, 'r') as f:
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
    return dic

def get_incar_dicts(files):
    dic_list = []
    for fname in files:
        if os.path.isdir(fname) and os.path.isfile(f"{fname}/INCAR"):
            fname = f"{fname}/INCAR"
        if os.path.isfile(fname):
            dic = get_incar_dict(fname)
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

def compare_incars(files, key, Ldiff=False):
    tup = get_incar_dicts(files)
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
                print_key(t, key)
            return 0

    d1keys = set(d1.keys())
    d2keys = set(d2.keys())
    shared = d1keys.intersection(d2keys)
    diff_1 = d1keys.difference(d2keys)
    diff_2 = d2keys.difference(d1keys)
    ### all the keys in common
    union  = d1keys.union(d2keys)
    ### sort the union in the order
    ordered_keys=[]
    for u in union:
        if not u in ordered_incar:
            print(f"Err: {u} is not listed in ordered_incar of mod_incar.py")
    for k in ordered_incar:
        if k in union:
            ordered_keys.append(k)
        
    s = ' '
    for key in ordered_keys:
        if key in shared:
            if d1[key] == d2[key]:
                diff="-"
                if Ldiff:
                    continue
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
   
def check_kw(ftype, kws):
    ### change kws into UPPERCASE
    kws = [ el.upper() for el in kws ]
    pwd = os.getcwd()
    files=[]
    if ftype == 'a' or ftype == 'd':
        allfiles = os.listdir(pwd)
        dirs = [ f for f in allfiles if os.path.isdir(pwd+'/'+f) ]
        dirs.sort()
        #print(f"{dirs}")
        for d in dirs:
            incar = pwd + '/'+d+'/INCAR'
            if os.path.isfile(incar):
                files.append(incar)
    if ftype == 'a' or ftype == 'f':
        allfiles = os.listdir(pwd)
        incars = [ f for f in allfiles if re.match('INCAR', f) ]
        incars.sort()
        files.extend(incars)
    print(f"{files}")
    i = 0
    #value_ini=[]
    for f in files:
        i += 1
        dic = get_incar_dict(f)
        #print(f"{dic} in {f}")
        if re.search('/', f):
            fs = re.split('/', f)
            print(f"{fs[-2]:20}:",end='')
        else:
            print(f"{f:20}:",end='')
        ### kws loop
        for j, kw in enumerate(kws):
            #if i == 1 :
                ### this makes error if dic[kw] does not exist
            #    value_ini.append(dic[kw])
            if kw in dic.keys():
                #if value_ini[j] == dic[kw]:
                print(f" {kw:10} {dic[kw]:>5} {'':10}", end='')
                #else:
                #    print(f" {kw:10} {dic[kw]:>5} {'diff':<10}", end='')
            else:
                print(f" {kw:10} {'None':>5}", end='')
        print("")
    return 0            

def main():

    parser = argparse.ArgumentParser(description="Compare two INCAR")
    parser.add_argument('-j', '--job', default='diff', choices=['diff','kw'], help="show INCAR difference or show keywords")
    parser.add_argument('-f', '--files', nargs='+',  help="fnames or symbol ['d','f','a']")
    parser.add_argument('-k', '--kws', nargs='*',  help="check the key and values ")
    parser.add_argument('-d', '--Ldiffer', action='store_true',  help="show only differences")
    parser.add_argument('-s', '--Lshow', action='store_true',  help="just show")
    args = parser.parse_args()

    file_choice=['d','f','a']
    if args.job == 'kw':
        if len(args.files) == 1 and args.files[0] in file_choice:
            if not args.kws:
                print(f"for -j kw -f [d,f,a] includes -kw kws")
            else:
                check_kw(args.files[0], args.kws) 
        else:
            print(f"for -j kw select {file_choice}")
    elif args.job == 'diff':
        compare_incars(args.files, args.kws, args.Ldiffer)

if __name__ == "__main__":
    main()
