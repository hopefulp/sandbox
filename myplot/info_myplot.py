#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_all, MyClass, dir_classify_n, whereami

my_mpl  = MyClass('my_mpl')
mplplot = MyClass('mplplot')
musage  = MyClass('musage')
amp     = MyClass('amp')

my_mpl.ini = "initialize mpl in ~/.config/matplotlib/matplotlibrc \
                \n\t\t:: check by ipython>>>matplotlib.matplotlib_fname()\
                "
mplplot.mplot_1f="myplot.py -v|-f values|files -j job -t title\
                \n\t\t:: -v y1 y2 y3 ... | -f f1 f2 f3 ...\
                \n\t\t:: -j for job qcmo|ai for xlabel, ylabel, title\
                \n\t\t:: -t, -xt, -yt overwrites xlabel, ylabel, title\
                \n\t\t:: -x for x-column -other options for title\
                \n\t\t--imports my_mplot2d for mplot_nvector\
                \n\t\t--imports plot_job for figure titles for jobs\
                \n\te.g.:(qcmo) myplot.py -v -y -0.8058   -0.7866   -0.7860   -1.0080   -1.2482   -1.2539 -j qcmo -t \"CO2 charges\"\
                \n\te.g.:(qcmo) myplot.py -f nbo-6f.dat -j qcmo -xt Model -ys -1 -yt \"Charge (e)\" -t \"CO2 charges\"\
                "
mplplot.mplot_f = mplplot.mplot_1f
mplplot.mplot_pd2f  =   "Draw 2d plot using pandas by reading file\
                        \n\tUsage:: mplot_pd2f.py fname options\
                        \n\tOptions::\
                        \n\t    -pd use pandas to ordering x-values\
                        \n\t\t-xs xspacing to reduce xticis by devide in the N list\
                        \n\t\t-xst ['evenly','numerically'(default)] in x-values\
                        "
mplplot.mplot_pdnf =    "extension of mplot_pd2f.py to include many files\
                        \n\tUsage:\
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
amp.md = " myplot.py md.ene -x -t MD-Ethylene -yt \"E(eV)\" -xt \"time (10fs)\" "
#classobj_dict={'MPL': my_mpl, 'MYPLOT': myplot, 'USAGE': musage}
classobj_dict={'MPL': my_mpl, 'MYPLOT': mplplot, 'USAGE': musage}

def classify(Lclassify, work, job):

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
        if not Lclassify:
            for f in sort_exe:
                print(f"    {f}")
        else:
            ### confer "ini_pypbs.py", MyClass.instances is a list of string as class variables
            for instance in MyClass.instances:
                ### globals() includes MyClass() instances as keys
                for gkey in globals().keys():
                    if gkey == instance.name:
                        break
                ### work == None without -w
                if work != instance.name:
                    ckeys = dir_classify_n(sort_exe, instance.name, globals()[gkey], Lwrite=1) # globals()[instance.name] is not working
                else:
                    ckeys = dir_classify_n(sort_exe, instance.name, globals()[gkey], Lwrite=0)
                    for ckey in ckeys:
                        print(f"    {ckey}.py\t:: {globals()[gkey].__dict__[ckey]}")
            print("  == not classified")
            for f in sort_exe:
                print(f"    {f}")

    if sort_mod:
        print("Module:: ")
        if not Lclassify:
            for f in sort_mod:
                print("    {}".format(f))
        else:
            for instance in MyClass.instances:
                for gkey in globals().keys():
                    if gkey == instance.name:
                        break
                if not work or work != instance.name:
                    ckeys = dir_classify_n(sort_mod, instance.name, globals()[gkey], Lwrite=1)
                else:
                    ckeys = dir_classify_n(sort_mod, instance.name, globals()[gkey], Lwrite=0)
                    for ckey in ckeys:
                        print(f"    {ckey}.py\t:: {globals()[gkey].__dict__[ckey]}")
            print("  == not classified ")
            for f in sort_mod:
                print(f"    {f}")


    ### print comment here
    if job in classobj_dict.keys():
        name_class = classobj_dict[job]
        for key in name_class.__dict__.keys():
            print(f" {job}   \t:: {name_class.__dict__[key]}")
    print("\nClass Instances:: ", end='')
    for instance in MyClass.instances:
        print(f"{instance.name}", end=' ')
    print("\n\t    -w for detail")

    return 0

def main():
    parser = argparse.ArgumentParser(description="display Usage for /mymplot  ")
    parser.add_argument('-c','--classify', action='store_false', help="classify ")
    parser.add_argument('-w','--work', help="explain not-file-related work")
    parser.add_argument('-cn','--classname', help="present class details ")
    #parser.add_argument('-js','--specify', choices=['qcmo','nbo','eda'], help="present class details ")
    parser.add_argument('-u','--usage', action='store_true', help="present main details")
    args = parser.parse_args()

    classify(args.classify, args.work, args.classname)

if __name__ == "__main__":
    main()
