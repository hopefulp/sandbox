#!/gpfs/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_files

def jobs(job):
    if job == None:
        print("List this directory = ")
    return 0 
    

def main():

    parser = argparse.ArgumentParser(description="display Usage for /mymplot  ")
    parser.add_argument('-j','--job',  help=" ")
    #parser.add_argument('-l','--list', action='store_true',  help="list directory files ")
    args = parser.parse_args()

    jobs(args.job)

if __name__ == "__main__":
    main()
