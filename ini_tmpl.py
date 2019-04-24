#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_all

source_dir = "/Research"
home = os.environ['HOME']
s_dir = home + source_dir

def jobs(job):
    
    if job == None:
        print("List this directory = ")
        mdir = s_dir
        exe, mod, dirs, d_link = dir_all(mdir)
        sort_exe = sorted(exe)
        sort_mod = sorted(mod)
        sort_dir = sorted(dirs)
        print("Directories:: ")
        for f in sort_dir:
            print(f"    {f}")
        print("Executable:: ")
        for f in sort_exe:
            print(f"    {f}")
        print("Module:: ")
        for f in sort_mod:
            print(f"    {f}")
        print("#Comment: ")
    

def main():
    parser = argparse.ArgumentParser(description="display Usage for /mymplot  ")
    parser.add_argument('-j','--job',  help=" ")
    #parser.add_argument('-l','--list', action='store_true',  help="list directory files ")
    args = parser.parse_args()

    jobs(args.job)

if __name__ == "__main__":
    main()
