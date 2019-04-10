#!/gpfs/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_files

amp_collection = {
    'file_conv':    'im2extxyz.py'  ,
    'amp_run':      'amp_ene.py'    ,
    'amp_valid':    'amp_validation.sh'}

models= {
    'ethylene': 'Ethylene.extxyz',
    'Diss_CHO': 'Diss_H2COH.extxyz' }

def fconvert_eg():
    print("    INF file to extxyz\n\t{} xxx.inf -a atom_list -y_bar [1,2,3]".format(amp_collection['file_conv']))
    print("\t: atom_list comes from 'molecules.py'")
    print("\t: -y_bar is options: 1 for energy only")
    return


def file_conversion():
    print("file conversion: {}".format(amp_collection['file_conv']))
    fconvert_eg()
    return 0


def run_amp():
    print("amp run        : {}".format(amp_collection['amp_run']))
    print("    For sample profile::\n\t{} {} profile".format(amp_collection['amp_run'],models['Diss_CHO']))
    print("    For job=tr(ain)::\n\t{} {} tr -hl 4 4 4 -el 0.001 -n 5".format(amp_collection['amp_run'],models['ethylene']))
    print("    For job=te(st)::\n\t{} {} te -hl 4 4 4 -el 0.001 -n 5".format(amp_collection['amp_run'],models['ethylene']))
    print("    For job=md::\n\t{} {} md".format(amp_collection['amp_run'],models['ethylene']))
    print("    For validation::\n\t{} -h\n\t{} {} '4 4 4' 0.001 +g | sh".format(amp_collection['amp_valid'],amp_collection['amp_valid'],models['ethylene']))
    return 0

def jobs(job):
    if job == None :
        print("List this directory = ")
        mdir = os.path.dirname(__file__)
        exe, mod = dir_files(mdir)
        sort_exe = sorted(exe)
        sort_mod = sorted(mod)
        print("Executable:: ")
        for f in sort_exe:
            print("    {}".format(f))
        print("Module:: ")
        for f in sort_mod:
            print("    {}".format(f))
        print("#Comment: try '-j amp'")
    elif job == 'amp':
        print("For AMP::")
        file_conversion()
        run_amp()

def main():

    parser = argparse.ArgumentParser(description="display Usage for /mymplot  ")
    parser.add_argument('-j','--job',  help=" ")
    #parser.add_argument('-l','--list', action='store_true',  help="list directory files ")
    #parser.add_argument('-ls','--list_detail', action='store_true',  help="list directory files ")
    args = parser.parse_args()

    jobs(args.job)

if __name__ == "__main__":
    main()
