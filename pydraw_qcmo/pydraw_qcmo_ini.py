#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_files

def jobs(job, li):
    if li:
        print("List this directory: ")
        mdir = os.path.dirname(__file__)
        exe, mod = dir_files(mdir)
        print("executable: {}".format(exe))
        print("module    : {}".format(mod))

def main():

    parser = argparse.ArgumentParser(description="display Usage for /mymplot  ")
    parser.add_argument('-j','--job',  help=" ")
    parser.add_argument('-l','--list', action='store_true',  help="list directory files ")
    args = parser.parse_args()

    jobs(args.job, args.list)

if __name__ == "__main__":
    main()
