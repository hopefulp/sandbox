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
                \n\t\t:: -j for job qcmo|ai for xlabel, ylabel, title\
                \n\t\t:: -t, -xt, -yt overwrites xlabel, ylabel, title\
                \n\t\t:: -x for x-column -other options for title\
                \n\t\t--imports my_mplot2d for mplot_nvector\
                \n\t\t--imports plot_job for figure titles for jobs\
                \n\te.g.:(qcmo) myplot.py -v -y -0.8058   -0.7866   -0.7860   -1.0080   -1.2482   -1.2539 -j qcmo -t \"CO2 charges\"\
                \n\te.g.:(qcmo) myplot.py -f nbo-6f.dat -j qcmo -xt Model -ys -1 -yt \"Charge (e)\" -t \"CO2 charges\"\
                "
musage.qcmo="(qcmo) myplot.py -v -y -0.8058   -0.7866   -0.7860   -1.0080   -1.2482   -1.2539 -j qcmo -t \"CO2 charges\"\
            \n\t\t\t  myplot.py -f nbo-6f.dat -j qcmo -xt Model -ys -1 -yt \"Charge (e)\" -t \"CO2 charges\" "
musage.eda="(eda) grep Polar *out | awk '{print $6}'\
            \n\t\t\t grep \"CT = DEL\" *out | awk '{print $11}' | tr '\\n' ' '\
            \n\t\t\t grep 'SCF Total' *out | awk '{print $11}' | tr '\\n' ' '\
            \n\t\t\t myplot.py -v -y -0.8058   -0.7866   -0.7860  -j eda -t 'CT Energy' -yt 'E (kcal/mol)' \
            \n\t\t\t myplot.py -f frozen_1.dat Polar.dat CTene.dat scf.dat -j eda -t EDA -yt 'E (kcal/mol)' -ys -1 -yl FRZ POL CT SCF-TOTAL\
            \n\t\t\t myplot.py -f chg-nbo.dat CTene.dat -ys -1 j- -yl 'NAO Charge of CO2 (e$^-$)' 'CT (kcal/mol)' -tx\
            \n\t\t\t myplot.py -f BE.dat scf.dat -ys -1 j- -t 'BE & SCF' -yt 'E (kcal/mol)' -yl BE SCF-TOTAL\
            \n\t\t\t myplot.py -f chg-nbo.dat BE.dat -ys -1 j- -yl 'NAO Charge of CO2 (e$^-$)' 'BE (kcal/mol)' -tx -c r darkcyan\
            \n\t\t\t myplot.py -f CTene.dat scf.dat -ys 'j-' 'j-' -yl 'CT (kcal/mol)' 'SCF (kcal/mol)' -tx -c red blue\
            "
amp     = MyClass()
amp.md = " myplot.py md.ene -x -t MD-Ethylene -yt \"E(eV)\" -xt \"time (10fs)\" "
classobj_dict={'MPL': my_mpl, 'MYPLOT': my_plot, 'USAGE': musage}
classobj_work={'AMP': amp}

def jobs(Lclass, job, work, Lusage):

    if Lusage:
        job = 'USAGE'
        name_class = classobj_dict[job]
        for key in name_class.__dict__.keys():
            print(f" {job}   \t:: {name_class.__dict__[key]}")
        return 0

    print("List this directory :: ")
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
    ### write for not-file-related work
    if not work:
        print(f"\n\t  -w {classobj_work.keys()} for not-file-related work")
    else:
        name_class = classobj_work[work]
        for key in name_class.__dict__.keys():
            print(f" {work}   \t:: {name_class.__dict__[key]}")

    return 0

def main():
    parser = argparse.ArgumentParser(description="display Usage for /mymplot  ")
    parser.add_argument('-c','--classify', action='store_true', help="classify ")
    parser.add_argument('-cn','--classname', help="present class details ")
    parser.add_argument('-w','--work', help="explain not-file-related work")
    #parser.add_argument('-js','--specify', choices=['qcmo','nbo','eda'], help="present class details ")
    parser.add_argument('-u','--usage', action='store_true', help="present main details")
    args = parser.parse_args()

    jobs(args.classify, args.classname, args.work, args.usage)

if __name__ == "__main__":
    main()
