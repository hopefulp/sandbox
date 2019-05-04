#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_all, MyClass

mo_exe = MyClass()
mo_exe.mplt_mo_dev="plot MO using Q-Chem output (job = sp)"



def jobs(job):
    
    if job == None:
        mdir = os.path.dirname(__file__)
        print(f"List directory of {mdir} ")
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
        print("#Comment: -j detail")
    elif re.search("de", job):
        print("Detail:: ")
        for key in mo_exe.__dict__.keys():
            print(f"    {key}\t:: {mo_exe.__dict__[key]}")
    

def main():
    parser = argparse.ArgumentParser(description="display Usage for $SB/py_qcmo  ")
    parser.add_argument('-j','--job',  help="classify ")
    #parser.add_argument('-l','--list', action='store_true',  help="list directory files ")
    args = parser.parse_args()

    jobs(args.job)

if __name__ == "__main__":
    main()
