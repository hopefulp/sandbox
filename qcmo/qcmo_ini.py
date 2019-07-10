#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_all, MyClass, dir_classify
from mplt_mo_ini import *

qcmo = MyClass()
mplt = MyClass()
nico2 = MyClass()
qchem = MyClass()
usage = MyClass()
qchem.version=  " <= v3.1 :: for MO plot of Jmol\
                \n\t\t\t           some key-word are not wokring\
                \n\t\t    v4.1    :: \
                \n\t\t\t          EDA runs for serial - No problem in multiplicity\
                \n\t\t    v5.1    ::\
                "


qcmo.mplot_mo="plot MO using Q-Chem output (job = sp)\n\tUsage    : mplot_mo.py -h\n\tcall     : mplt_ab_draw for drawing\
                    \n\t\t Block is modulized in qcout_mod.py\n\t \
                    \n\te.g. (level)   : mplot_mo.py -f 1-PP-A.out 1-PP.out 1-PP-B.out \
                    \n\te.g. (link 3file)   : mplot_mo.py -f 1-PP-A.out 1-PP.out 1-PP-B.out \
                    \n\t\t\t\t-a 'Ni' 'Ni C1 O1 O2' 'C O1 O2' -t ONE SUB ALL -l [-lf m1-3f.dat] \
                    \n\te.g. (link 5file)   : mplot_mo.py -f 1-PP-A.out 1-PP-A.out 1-PP.out 1-PP-B.out CO2.out \
                    \n\t\t\t\t-a 'Ni' 'Ni' 'Ni C1 O1 O2' 'C O1 O2' 'C O1 O2' -t ONE ONE SUB ALL ALL -l -lf m1-5f.dat  'w. link file'\
                    \n\tParameters::\
                    \n\t\t-l    : draw link, this makes 'link_id.dat', then save to Model1-3f.dat to modify link\
                    \n\t\t-lf   : use link index file, which has 'link1\\nlink2\\nlin...' etc\
                    \n\t\t-a    : atom types to select MO based on the selected atoms_indices -a 'Ni C1 O1 O2' etc\
                    \n\t\t-t --type: how to choose MO based on motype \
                    \n\t\t\tONE if moc-line has one atom in the atom list given by  \
                    \n\t\t\tSEL if moc-line has any of atom \
                    \n\t\t\tALL if moc-line has all the atoms \
                    "
qcmo.mplt_mo_ini="initial condition for plot: V_print fore verbose "
qcmo.qcout_mod="modularlized Blocks of MO Energies, MO Coefficients, NBO charges \
                \n\t\t\tgeometric average of coefficients of the same base is done at imo_basis_dic()"
qcmo.mo_level_link="functions::\n\t\t\tget_link(): called by main(); obtain link_id between two files; activated by -l [ltypes] by main\
                    \n\t\t\ttweak HOMO LUMO and link"    
mplt.mplt_qcdraw="called by mplot_mo.py through mplt_ab_draw()\
                \n\t\tset matplotlib.rcParams; cal x-axis; draw by mplt_ab_draw\n\t\t\tHelp:: python $SB/qcmo/mplt_qcdraw.py\
                \n\t\tmplt_ab_draw():: calls\
                \n\t\t\tmplt_level():: calls\
                \n\t\t\t\tmplot_arow_homo():: draw arrow for occupied levels\
                \n\t\t\tXrange_nf_fixed_x_length():: use decimal.Decimal for x[]"

nico2.Ni_CO2red="Usage:: Ni_CO2red.py -j {moc,nbo} -m {1,2,3,4,5} -nf {1, 3, 5} -l -lf -bar -r\
                \n\t\t-j : job of MOC and NBO\
                \n\t\t  -js, --job_level 2 [default=1] for more extension\
                \n\t\t-m : model 1:1-PP, 2:2-PPP, 3:3-M3-PNP, 4:4-PNP, 5:5-FeNi, 6:6-CC-FeNi\
                \n\t\t-nf : number of qcout files 3: m-A m m-B, 5: m-A-relax 3-files CO2-relax\
                \n\t\t-l : make link by calculation\
                \n\t\t-lf : use link_id.dat\
                \n\t\t-bar : use m-A-relax instead of frozen m-A\
                \n\t\t-r : ask whether run or copy the command\
                \n\t\te.g.: Ni_CO2red.py -m 4 -nf 5 "
nico2.usage="\tNi_CO2.py -j moc -m 4 -nf 5 -l -lf\
            \n\tNi_CO2red.py -j nbo -js 2 "
usage.mlot_mo_mod = qcmo.mplot_mo

classobj_dict={'QCMO':qcmo, 'MPLT':mplt, 'NiCO2':nico2, 'QChem':qchem, 'Usage':usage}

def jobs(Lclass,job, Lusage, model):
    if Lusage:
        job = 'Usage'
        name_class = classobj_dict[job]
        for key in name_class.__dict__.keys():
            print(f" {job}   \t:: {name_class.__dict__[key]}")
        return 0

    mdir = os.path.dirname(__file__)
    print(f"List directory of {mdir} ")
    exe, mod, dirs, d_link = dir_all(mdir)
    sort_exe = sorted(exe)
    sort_mod = sorted(mod)
    sort_dir = sorted(dirs)
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
    if job == "QChem":
        name_class = classobj_dict[job]
        for key in name_class.__dict__.keys():
            print(f" {job}   \t:: {name_class.__dict__[key]}")
            
    print(f"#Comment: -c    for classification'\
            \n\t  -u for usage: equal to -j USAGE\
            \n\t  -j {classobj_dict.keys()} for detail ")

def main():
    parser = argparse.ArgumentParser(description="display Usage for $SB/py_qcmo  ")
    parser.add_argument('-c','--classify', action='store_true', help="classify ")
    parser.add_argument('-j','--job', help="present class details ")
    parser.add_argument('-u','--usage', action='store_true', help="present main details")
    parser.add_argument('-m','--model', type=int, help="present main details")
    args = parser.parse_args()

    jobs(args.classify,args.job,args.usage,args.model)

if __name__ == "__main__":
    main()
