#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_all, MyClass

mo_exe = MyClass()
mo_exe.mplt_mo_dev="plot MO using Q-Chem output (job = sp)\n\tUsage    - mplt_mo_dev.py -h"

mo_mod = MyClass()

classobj_dict={'EXE':mo_exe, "MOD":mo_mod}

def jobs(job,spec):
    
    if job == None or re.search("cl", job):
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
        if job == None:
            for f in sort_mod:
                print(f"    {f}")
        else:
            luse=[]
            for f in sort_mod:
                f1=re.split(".",f)[0]
                if f1 in mo_mod.__dict__.keys():
                    luse.append(f)
                    continue
                print("    {}".format(f))
            ### classify modules used
            print("  {:<10}::".format('MOD used'))
            for f in luse:
                print(f"{f}")
                    
        print("#Comment: -j classify   : for more classification")
        print("\t:  -s Class\t: for specification for the class")
        if not spec == None:
            print(f"Detail for {spec}:: ")
            name_class = classobj_dict[spec]
            for key in name_class.__dict__.keys():
                print(f"    {key}\t:: {name_class.__dict__[key]}")
    

def main():
    parser = argparse.ArgumentParser(description="display Usage for $SB/py_qcmo  ")
    parser.add_argument('-j','--job',  help="classify ")
    parser.add_argument('-s','--specification',  help="classify ")
    args = parser.parse_args()

    jobs(args.job,args.specification)

if __name__ == "__main__":
    main()
