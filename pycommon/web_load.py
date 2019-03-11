#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import yes_or_no

def run(dir1,dir2,files,L_run):   
    load_files=[]
    home = os.getenv('HOME')
    if not dir1:
        dir1 = "/shared/share_win/web_load/"
    if not dir2:
        dir2 = home + "/public_html/research/"

    if not files:
        load_files = os.listdir(dir1)
    else:
        load_files = files

    for f in load_files:
        fname =  dir1 + f    
        com = "cp %s %s" % (fname, dir2)
        if L_run or yes_or_no(com):
            os.system(com)

    return 0

def main():
    parser = argparse.ArgumentParser(description='web-load files')
    parser.add_argument('-d1', '--ini', help='directory file sources')
    parser.add_argument('-d2', '--destiny', help='directory for web loading')
    parser.add_argument('-f', '--files', nargs='*', help='specify files to be loaded')
    parser.add_argument('-r', '--run', action='store_true', help='run without asking')
    args = parser.parse_args()
    
    run(args.ini, args.destiny,args.files, args.run)
    return 0

if __name__=='__main__':
	main()	


