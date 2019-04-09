#!/home/joonho/anaconda3/bin/python

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
        for f in exe:
            print("    {}".format(f))
        print("Module:: ")
        for f in mod:
            print("    {}".format(f))
        print("#Comment: try '-j myplot'")
    elif job == 'myplot':
        print("myplot.py file[files] -x for x-column -other options for title")
        print("Plot MD: \n    myplot.py md.ene -x -t MD-Ethylene -yt \"E(eV)\" -xt \"time (10fs)\" ")

def main():

    parser = argparse.ArgumentParser(description="display Usage for /mymplot  ")
    parser.add_argument('-j','--job',  help=" ")
    #parser.add_argument('-l','--list', action='store_true',  help="list directory files ")
    args = parser.parse_args()

    jobs(args.job)

if __name__ == "__main__":
    main()
