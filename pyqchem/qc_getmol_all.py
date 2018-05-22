#!/home/jackjack5/epd/bin/python

import argparse
import os
import re
import subprocess as subp
from common import *

#sh_cmd="/qcfs/joonho/bin/qgeo_mo.pl"

def get_mol(ftype, job_type, prefix):
    pwd = os.getcwd()
    suff = []
    suff.append(ftype)     # get_files_suffix receives list of suffix
    ### obtain files of suffix |in|out| from directory 
    f_list = get_files_suffix(suff, pwd)
    f_ts = []
    f_opt = []
    #if not job_type:
    #    job_type == "opt"
    #    
    #if job_type == "opt":
    #    n = 1
    #elif job_type == "ts":
    #    n = 3
    
    ### divide files into |opt|ts|
    for file in f_list:
        if re.search("ts", file, re.IGNORECASE):
            f_ts.append(file)
        else:
            f_opt.append(file)
    # in case not all the files
    #if prefix:
    #    if not re.match(prefix, file):
    #        continue
    #    else:
    #        print "match ", file
    #if not re.search("out$", file):
    if not job_type or job_type == "opt":
        for file in f_opt:
            cmd = "qgeo_mo.pl %s 1 " % file 
            print cmd
            os.system(cmd)
    if not job_type or job_type == "ts":
        for file in f_ts:
            cmd = "qgeo_mo.pl %s 3 " % file
            print cmd
            os.system(cmd)
    # print cmd
    #if fout_type == "ts":
    #    subp.call([sh_cmd, file,  "3"])
    #else:
    #    subp.call([sh_cmd, file])



def main():
    parser = argparse.ArgumentParser(description='obtain mol file from Q-Chem file: opt, ts')
    parser.add_argument("-f", "--file_type", type=str, default="out", choices=["out", "in"], help='default file type is .out')
    parser.add_argument("-p", "--prefix", type=str, help='prefix for filetype')
    parser.add_argument("-j", "--jobtype", choices=["opt","ts","sp"], help='qchem job-type')
    args = parser.parse_args()
    get_mol(args.file_type, args.jobtype, args.prefix)


if __name__ == "__main__":
    main()
