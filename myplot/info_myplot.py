#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_all, MyClass, dir_classify_n, whereami
#from comment_subj import vasp

my_mpl  = MyClass('my_mpl')
mplplot = MyClass('mplplot')
usage   = MyClass('usage')
table   = MyClass('table')
test    = MyClass('test')

my_mpl.ini = "initialize mpl in ~/.config/matplotlib/matplotlibrc \
                \n\t\t::usage.amp check by ipython>>>matplotlib.matplotlib_fname()\
                "
mplplot.mplot_f="myplot.py -v|-f values|files -j job -t title\
                \n\t\t:: -v y1 y2 y3 ... | -f f1 f2 f3 ...\
                \n\t\t:: -j for job qcmo|ai for xlabel, ylabel, title\
                \n\t\t:: -t, -xt, -yt overwrites xlabel, ylabel, title\
                \n\t\t:: -x for x-column -other options for title\
                \n\t\t--imports my_mplot2d for mplot_nvector\
                \n\t\t--imports plot_job for figure titles for jobs\
                \n\te.g.:(qcmo) myplot.py -v -y -0.8058   -0.7866   -0.7860   -1.0080   -1.2482   -1.2539 -j qcmo -t \"CO2 charges\"\
                \n\te.g.:(qcmo) myplot.py -f nbo-6f.dat -j qcmo -xt Model -ys -1 -yt \"Charge (e)\" -t \"CO2 charges\"\
                "
mplplot.mplot_pdnf =    "plot multiple files using pandas for ordering\
                        \n\t    any number of files\
                        \n\t    pandas is default\
                        \n\t    filename becomes columns Name\
                        \n\t    mplot_pd2f option is available\
                        \n\tUsage::\
                        \n\t    mplot_pdnf.py DataAnal/tr300.dat DataAnal/tr500.dat DataAnal/te_tr300.dat\
                        \n\tOptions::\
                        \n\t    -pd is default to use pandas to ordering x-values\
                        \n\t\t-xs xspacing to reduce xticis by devide in the N list\
                        \n\t\t-xst ['evenly','numerically'(default)] in x-values\
                        \n\t(mplot_pdnf2) Draw 2d plot using pandas by reading file\
                        \n\t    will be deprecated to use mplot_pdnf.py\
                        \n\tUsage:: mplot_pd2f.py 1-2filenames options\
                        "
mplplot.myplot2D    =   "several kinds of 2d plot method\
                        \n\tdef mplot_twinx:\
                        \n\tdef mplot_nvector:\
                        \n\t    called by:\
                        \n\t\tmplot_f.py\
                        \n\t\tmplot_table.py\
                        \n\t\tplot_level.py\
                        \n\t\trun_plot_qcmo.py\
                        \n\tdef mplot_vector_one:\
                        \n\tdef mplot_vector_two:\
                        \n\tdef draw_histogram:\
                        \n\tdef barplot2:\
                        \n\tdef barplot_y:\
                        "
table.mplot_gibbs   =   "plot csv: written by David Park\
                        \n\tHow to Use\
                        \n\t    table w. white space without empty (use nan in empty space) works with 'mv a.dat a.csv'\
                        \n\t    save excel sheet to csv format\
                        \n\tUsage\
                        \n\t    $mplot_gibbs.py MXene-4level.csv -l 'G(U=0)' 'G(\$U_{Dc}$=1.03)' 'G(\$U_{Eq}$=2.79)' 'G(\$U_{Ch}$=4.79)' -c k b g r\
                        \n\t    $mplot_gibbs.py MXene-4level.csv -l 'G(U=0)' 'G($U_{Dc}$=1.03)' 'G($U_{Eq}$=2.79)' 'G($U_{Ch}$=4.79)' -xl 'O$_2$ Reduction step' -c k b g r\
                        \n\t    $mplot_gibbs.py apcc.csv -c r r b b k\
                        \n\tText\
                        \n\t    hypertext: $_{}$ $^{}$ (in chi) or \\$_{}$ \\$^{}$ (in Pt)\
                        "
### this should follow the definition of table.mplot_gibbs
mplplot.mplot_gibbs = table.mplot_gibbs                        
table.mplot_table   =   "Plot table: CSV from MS Excel\
                        \n\t    import plot_level.mplot_level\
                        "
table.plot_level    =   "modules for table plot\
                        \n\t    imported in mplot_table\
                        \n\t    modified from mplot_f.py\
                        "
usage.qcmo      = "myplot.py -v -y -0.8058   -0.7866   -0.7860   -1.0080   -1.2482   -1.2539 -j qcmo -t 'CO2 charges'\
                \n\t\t\t  myplot.py -f nbo-6f.dat -j qcmo -xt Model -ys -1 -yt 'Charge (e)' -t 'CO2 charges'\
                "
usage.eda       = "grep Polar *out | awk '{print $6}'\
                \n\t\t\t grep \"CT = DEL\" *out | awk '{print $11}' | tr '\\n' ' '\
                \n\t\t\t grep 'SCF Total' *out | awk '{print $11}' | tr '\\n' ' '\
                \n\t\t\t myplot.py -v -y -0.8058   -0.7866   -0.7860  -j eda -t 'CT Energy' -yt 'E (kcal/mol)' \
                \n\t\t\t myplot.py -f frozen_1.dat Polar.dat CTene.dat scf.dat -j eda -t EDA -yt 'E (kcal/mol)' -ys -1 -yl FRZ POL CT SCF-TOTAL\
                \n\t\t\t myplot.py -f chg-nbo.dat CTene.dat -ys -1 j- -yl 'NAO Charge of CO2 (e$^-$)' 'CT (kcal/mol)' -tx\
                \n\t\t\t myplot.py -f BE.dat scf.dat -ys -1 j- -t 'BE & SCF' -yt 'E (kcal/mol)' -yl BE SCF-TOTAL\
                \n\t\t\t myplot.py -f chg-nbo.dat BE.dat -ys -1 j- -yl 'NAO Charge of CO2 (e$^-$)' 'BE (kcal/mol)' -tx -c r darkcyan\
                \n\t\t\t myplot.py -f CTene.dat scf.dat -ys 'j-' 'j-' -yl 'CT (kcal/mol)' 'SCF (kcal/mol)' -tx -c red blue\
                "
usage.sno2      = "pyvasp/\
                \n\tdoslm.py -z 3.69\
                "
usage.h2        = "H2 on Pt-C60-x\
                \n\tf{vasp.scripts.zpe}\
                "
usage.amp       = " myplot.py md.ene -x -t MD-Ethylene -yt \"E(eV)\" -xt \"time (10fs)\" \
                \n\t\t\tmplot_f.py -v 1 2 3 4\
                "
usage.test      = "\t(ascii file) Test for mpl plot\
                \n\t    Usage:: mplot_f.py -v 1 2 3 4\
                \n\t\t    mplot_f.py -f t.dat\
                \n\t(csv file) Test for table plot\
                \n\t    Usage:: mplot_table.py GC54Pt-abc.csv -l -t 'H2 diffusion on C60-x'\
                "
#classobj_dict={'MPL': my_mpl, 'MYPLOT': myplot, 'USAGE': musage}
classobj_dict={'MPL': my_mpl, 'MYPLOT': mplplot, 'USAGE': usage}

def classify(Lclassify, work, classname, job):

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
    if work == 'usage':
        print("==== USAGE ====")
        for key in usage.__dict__.keys():
            print(f"{key.upper()}: {usage.__dict__[key]}")
    return 0

def main():
    parser = argparse.ArgumentParser(description="display Usage for /mymplot  ")
    parser.add_argument('-c','--classify', action='store_false', help="classify ")
    parser.add_argument('-w','--work', help="explain not-file-related work")
    parser.add_argument('-cn','--classname', help="present class details ")
    parser.add_argument('-j','--job', choices=['qcmo','nbo','eda'], help="present class details ")
    args = parser.parse_args()

    classify(args.classify, args.work, args.classname, args.job)

if __name__ == "__main__":
    main()
