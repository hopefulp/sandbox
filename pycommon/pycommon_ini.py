#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_all, MyClass, dir_classify

mod_dir = MyClass()
mod_dir.common="module with commonly used functions \"import common\""


exe_dir= MyClass()
exe_dir.dir_clean_p2="clean dir by -prefix -suffix -middle match -e excluded -ef 'exclude these files' -work {qchem,ai} -j rm[mv] -jd new_dir\n\t\tUsage   :: dir_clean_p2.py -s out -ef 6-CC-NiFe-A-relax.out 5-FePNP-CO2.out -j mv -jd j631gs_v3.2"
exe_dir.dir_clean_rec_p2="clean dir recursively by -p -s -m"
exe_dir.dir_reset="reset dir as initial state by job: -j ai"
exe_dir.dir_run="scan dir and run the same command for all the files such as\n\t\tqcout_mol_in.pl"

classobj_dict={'EDIR':exe_dir, "MDIR":mod_dir}

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
        if job == None:
            for f in sort_exe:
                print(f"    {f}")
        else:
            dir_classify(sort_exe, 'EDIR', classobj_dict)
        print("Module:: ")
        if job == None:
            for f in sort_mod:
                print(f"    {f}")
        else:
            dir_classify(sort_mod, 'MDIR', classobj_dict)
                    
        print("#Comment: -j classify   : for more classification")
        print("\t:  -s Class\t: for specification for the class")
        print("\t: *_p2.py  \t: python 2.7 ")
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
