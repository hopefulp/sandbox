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
        exe_s = sorted(exe)
        mod_s = sorted(mod)
        for f in exe_s:
            print("    {}".format(f))
        print("Module:: ")
        for f in mod_s:
            print("    {}".format(f))
        print("#Comment: \n    sge_home.sh fname.py\n\tmakes sge_fname.py with \"#!/gpfs/home/...\"")

    

def main():

    parser = argparse.ArgumentParser(description="display Usage for /mymplot  ")
    parser.add_argument('-j','--job',  help=" ")
    #parser.add_argument('-l','--list', action='store_true',  help="list directory files ")
    args = parser.parse_args()

    jobs(args.job)

if __name__ == "__main__":
    main()
