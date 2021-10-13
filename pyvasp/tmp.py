#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_files

def compare_incar(files):
    dic_list = []
    for fname in files:
        dic={}
        if os.path.isfile(fname):
            with open(fname, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    strline=line.strip()
                    if re.match('\w',strline):
                        if re.search('=', strline):
                            #prog = re.compile('[=\s]')
                            prog='[=\s]+'
                            lst = re.split(prog, strline)
                            print(f'{lst[0]:^10} = {lst[1]:>10}')
                            dic[lst[0]] = lst[1]
        dic_list.append(dic)                            
    return 0 
    

def main():

    parser = argparse.ArgumentParser(description="Compare two INCAR")
    parser.add_argument('files', nargs='+',  help="the 1st INCAR file ")
    #parser.add_argument('f2',  help="the 2nd INCAR file ")
    args = parser.parse_args()

    compare_incar(args.files)

if __name__ == "__main__":
    main()
