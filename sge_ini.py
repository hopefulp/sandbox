#!/gpfs/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_files

def jobs(job):
    if job == None:
        print("List this directory = ")
        mdir = os.path.dirname(__file__)
        exe, mod = dir_files(mdir)
        print("Executable:: ")
        sort_exe = sorted(exe)
        sort_mod = sorted(mod)
        for f in sort_exe:
            print("    {}".format(f))
        print("Module:: ")
        for f in sort_mod:
            print("    {}".format(f))
        print("#Comment: Here is SGE server in mlet")
        print("    SGE::(start)\n\tsge_pypbs_ini.py")
        
    

def main():

    parser = argparse.ArgumentParser(description="display Usage for /mymplot  ")
    parser.add_argument('-j','--job',  help=" ")
    #parser.add_argument('-l','--list', action='store_true',  help="list directory files ")
    args = parser.parse_args()

    jobs(args.job)

if __name__ == "__main__":
    main()
