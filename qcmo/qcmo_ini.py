#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_all, MyClass

mo_exe = MyClass()
mo_exe.mplt_mo_dev="plot MO using Q-Chem output (job = sp)\n\tUsage    : mplt_mo_dev.py -h\n\tcall     : mplt_ab_draw for drawing"

mo_mod = MyClass()
mo_mod.mplt_qcdraw="set matplotlib.rcParams; cal x-axis; draw by mplt_ab_draw\n\tHelp:: python $SB/qcmo/mplt_qcdraw.py"

classobj_dict={'EXE':mo_exe, "MOD":mo_mod}

def f_classify(lsorted, classobj_dict_key):
    #print(classobj_dict_key)
    c_obj = classobj_dict[classobj_dict_key]
    luse=[]
    for f in lsorted:
        f1=re.split("\.",f)[0]
        if f1 in c_obj.__dict__.keys():
            luse.append(f)
            continue
        print(f"    {f}")
    ### classify modules used
    print("  {:<10}::".format(classobj_dict_key+" used"))
    for f in luse:
        print(f"    {f}")


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
            f_classify(sort_exe, 'EXE')
            
        print("Module:: ")
        if job == None:
            for f in sort_mod:
                print(f"    {f}")
        else:
            f_classify(sort_mod, 'MOD')
        print("#Comment: -j classify   : for more classification")
        print("\t:  -s Class\t: for specification for the class")
        if not spec == None:
            print(f"Detail for {spec}:: ")
            name_class = classobj_dict[spec]
            for key in name_class.__dict__.keys():
                print(f"    {key}.py\t:: {name_class.__dict__[key]}")
    

def main():
    parser = argparse.ArgumentParser(description="display Usage for $SB/py_qcmo  ")
    parser.add_argument('-j','--job',  help="classify ")
    parser.add_argument('-s','--specification',  help="classify ")
    args = parser.parse_args()

    jobs(args.job,args.specification)

if __name__ == "__main__":
    main()
