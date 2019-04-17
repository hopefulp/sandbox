#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_files

amp_collection = {
    'file_conv':    'im2extxyz.py'  ,
    'amp_run':      'amp_ene.py'    ,
    'amp_valid':    'amp_validation.sh',
    'amp_scan':     'amp_loop.sh'}

models= {
    'ethylene': 'Ethylene.extxyz',
    'Diss_CHO': 'Diss_H2COH.extxyz',
    'water'   : 'water128.extxyz'}

def fconvert_eg():
    print("    INF file to extxyz\n\t{} xxx.inf -a atom_list -y_bar [1,2,3]".format(amp_collection['file_conv']))
    print("\t: atom_list comes from 'molecules.py'")
    print("\t: -y_bar is options: 1 for energy only")
    return


def file_conversion():
    print("file conversion: {}".format(amp_collection['file_conv']))
    fconvert_eg()
    return 0


def run_amp(fname,HL,elimit,nc,Lgraph):
    hl = " ".join(HL)
    print("amp run        : {}".format(amp_collection['amp_run']))
    print("    For sample profile::\n\t{} {} profile".format(amp_collection['amp_run'],models['Diss_CHO']))
    print(f"    For job=tr(ain)::\n\t{amp_collection['amp_run']} {fname} tr -hl {hl} -el {elimit} -n 5 -nc {nc}")
    if Lgraph:
        print(f"    For job=te(st)::\n\t{amp_collection['amp_run']} {fname} te -hl {hl} -el {elimit} -n 5 -nc {nc} +g")
    else:
        print(f"    For job=te(st)::\n\t{amp_collection['amp_run']} {fname} te -hl {hl} -el {elimit} -n 5 -nc {nc} -g")
    print("    For job=md::\n\t{} {} md".format(amp_collection['amp_run'],models['ethylene']))
    print("    For validation::\n\t{} -h\n\t{} {} '4 4 4' 0.001 +g | sh".format(amp_collection['amp_valid'],amp_collection['amp_valid'],models['ethylene']))
    print("    For validation scan::")
    print("                    [scan|not] job-type[val|tr] fname n_core")
    print("\t{} -h\n\t{} scan val {} 6 ".format(amp_collection['amp_scan'],amp_collection['amp_scan'],models['water']))
    return 0

def jobs(job,fname,HL, elimit, nc, Lgraph):
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
        run_amp(fname,HL, elimit, nc, Lgraph)

def main():

    parser = argparse.ArgumentParser(description="display Usage for ~/py_ai")
    parser.add_argument('-j','--job',  help="[val,train,test] ")
    parser.add_argument('-f','--file',  help="input energy data file ")
    parser.add_argument('-hl','--hidden_layer',nargs='*', default=['4','4','4'], help="list of number of Hidden Layer")
    parser.add_argument('-el','--energy_limit',default=0.001, type=float,  help="energy_limit for training")
    parser.add_argument('-nc','--ncore',  help="number of parallel process")
    parser.add_argument('-g','--graph', action='store_true',  help="draw graph or not")
    args = parser.parse_args()

    jobs(args.job, args.file, args.hidden_layer, args.energy_limit,args.ncore, args.graph )

if __name__ == "__main__":
    main()
