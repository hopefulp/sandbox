#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_all, MyClass, dir_classify

my_mpl = MyClass()
my_plot = MyClass()
musage = MyClass()
my_mpl.ini = "initialize mpl in ~/.config/matplotlib/matplotlibrc \
                \n\t\t:: check by ipython>>>matplotlib.matplotlib_fname()\
                "
my_plot.myplot="myplot.py -v|-f values|files -j job -t title\
                \n\t\t:: -v y1 y2 y3 ... | -f f1 f2 f3 ...\
                \n\t     (md  ) myplot.py md.ene -x -t MD-Ethylene -yt \"E(eV)\" -xt \"time (10fs)\" \
                "
musage.qcmo="(qcmo) myplot.py -v -y -0.8058   -0.7866   -0.7860   -1.0080   -1.2482   -1.2539 -j qcmo -t \"CO2 charges\"\
            \n\t\t\t  myplot.py -f nbo-6f.dat -j qcmo -xt Model -ys -1 -yt \"Charge (e)\" -t \"CO2 charges\" "
musage.eda="(eda) grep Polar *out | awk '{print $6}'\
            \n\t\t\t myplot.py -f CTene.dat scf.dat -ys 'j-' 'j-' -yl 'CT (kcal/mol)' 'SCF (kcal/mol)' -tx -c red blue\
            "


classobj_dict={'MPL': my_mpl, 'MYPLOT': my_plot, 'USAGE': musage}

def jobs(Lclass, job, Lusage):

    if Lusage:
        job = 'USAGE'
        name_class = classobj_dict[job]
        for key in name_class.__dict__.keys():
            print(f" {job}   \t:: {name_class.__dict__[key]}")
        return 0

    print("List this directory = ")
    mdir = os.path.dirname(__file__)        # __file__: this file location
    exe, mod, dirs, d_link = dir_all(mdir)
    sort_exe = sorted(exe)
    sort_mod = sorted(mod)
    sort_dir = sorted(dirs)
    #print(f"{exe} {mod}")
    if sort_dir:
        print("Directories:: ")
        for f in sort_dir:
            print(f"    {f}")

    if sort_exe:
        print("Executable:: ")
        if not Lclass:
            for f in sort_exe:
                print(f"    {f}")
        else:
            for key in classobj_dict.keys():
                fkey = dir_classify(sort_exe, key, classobj_dict)
                if job:
                    name_class = classobj_dict[job]
                    for key in name_class.__dict__.keys():
                        if key in fkey:
                            print(f"    {key}.py\t:: {name_class.__dict__[key]}")
            print("  == remainder ")
            for f in sort_exe:
                print(f"    {f}")
    if sort_mod:       
        print("Module:: ")
        if not Lclass:
            for f in sort_mod:
                print(f"    {f}")
        else:
            for key in classobj_dict.keys():
                fkey = dir_classify(sort_mod, key, classobj_dict)
                if job:
                    name_class = classobj_dict[job]
                    for key in name_class.__dict__.keys():
                        if key in fkey:
                            print(f"    {key}.py\t:: {name_class.__dict__[key]}")
            print("  == remainder ")
            for f in sort_mod:
                print(f"    {f}")
    if job == "MPL" or job == "USAGE":
        name_class = classobj_dict[job]
        for key in name_class.__dict__.keys():
            print(f" {job}   \t:: {name_class.__dict__[key]}")
    print(f"#Comment: -c    for classification'\
            \n\t  -u for usage: equal to -j USAGE\
            \n\t  -j {classobj_dict.keys()} for detail\
            ")


    return 0

def main():
    parser = argparse.ArgumentParser(description="display Usage for /mymplot  ")
    parser.add_argument('-c','--classify', action='store_true', help="classify ")
    parser.add_argument('-j','--job', help="present class details ")
    #parser.add_argument('-js','--specify', choices=['qcmo','nbo','eda'], help="present class details ")
    parser.add_argument('-u','--usage', action='store_true', help="present main details")
    args = parser.parse_args()

    jobs(args.classify, args.job, args.usage)

if __name__ == "__main__":
    main()
