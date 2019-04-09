#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_files

def jobs(job, li, li_detail):
    if li or li_detail:
        print("List this directory:: ")
        mdir = os.path.dirname(__file__)
        exe, mod = dir_files(mdir)
        if not li_detail:
            print("executable: {}".format(exe))
            print("module    : {}".format(mod))
        else:
            print("Executable:: ")
            for f in exe:
                print("    {}".format(f))
            print("Module:: ")
            for f in mod:
                print("    {}".format(f))


def main():

    parser = argparse.ArgumentParser(description="display Usage for /mymplot  ")
    parser.add_argument('-j','--job',  help=" ")
    parser.add_argument('-l','--list', action='store_true',  help="list directory files ")
    parser.add_argument('-ls','--list_detail', action='store_true',  help="list directory files ")
    args = parser.parse_args()

    jobs(args.job, args.list,args.list_detail )

if __name__ == "__main__":
    main()
