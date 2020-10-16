#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_all, MyClass, dir_classify_n, whereami

comment = MyClass('comment')
mod     = MyClass('mod')
dirjob  = MyClass('dirjob')
jobfile = MyClass('jobfile')
convert = MyClass('convert')
command = MyClass('command')

comment.howto       =   "=== COMMENT ===\
                        \n\treads 'comment_sys.py' or 'comment_subj.py'\
                        "
comment.showall     =   "runs 'howto.py'\
                        \n\tshows comment for system info\
                        "
comment.comment_sys =   "\tINFO: SYSTEM\
                        \n\tshowall.py\
                        "
comment.comment_subj=   "\tINFO: SCIENCE JOB \
                        \n\tshowall.py -s\
                        "
comment.ini_pycommon=   "\texplanation of files in $SB/pycommon"

mod.common="module with commonly used functions \"import common\""


dirjob.clean    =   "=== DIRECTORY JOB ==\
                        \n\t\tclean directory\
                        \n\t\t    by -prefix -suffix -middle match -e excluded -ef 'exclude these files' -work {qchem,ai} -j rm[mv] -jd new_dir\
                        \n\t\tUsage::\
                        \n\t\t    dir_clean_p2.py -s out -ef 6-CC-NiFe-A-relax.out 5-FePNP-CO2.out -j mv -jd j631gs_v3.2\
                        "
dirjob.dir_clean_r_py2= "clean dir recursively\
                        \n\t\tdir_clean_r_py2.py -p -s -m\
                        "
dirjob.dir_cli      =   "basic command line interface for all files in directory\
                        \n\t\tmodify script for all the files/selected files\
                        \n\t\tdir_cli.sh [0:vmake 1:incar 2:qsub 3:cp 4:rm 5:chmod\
                        "
dirjob.dir1_cli     =   "simple command inside directory\
                        \n\t\tdir1_cli.sh gitpush\
                        "
dirjob.diramp       =   "Run multiple job in amp by scanning a value in bash\
                        \n\t\tUsage::\
                        \n\t\t    $diramp.sh job other-args\
                        \n\t\tOptions::\
                        \n\t\t    fp: calculate fingerprints\
                        \n\t\t\tamp_wrapper.py -js qsub -j tr -qn NN9p\$num -dl \$num  \$(expr \$num + 360) &\
                        \n\t\t    wrapper: run amp_wrapper.py in multiple dires\
                        \n\t\t\tamp_wrapper.py -js qsub -qn \$dirname -k \$n & \
                        \n\t\t    te: run amp_wrapper.py -js qsub -j te\
                        "
dirjob.dir_fname    =   "jobs for ls, mvdir, rm, rename, cp\
                        \n\t\tOptions::\
                        \n\t\t    -[p|s|m] for matching type\
                        \n\t\t    -rp for replacement of matching part\
                        \n\t\t    -id to include dir in scanning dir\
                        \n\t\tUsage::\
                        \n\t\t    dir_fname.py rename -p G4 -rp G2 -id              ! to rename directories\
                        \n\t\t    dir_fname.py rm -m '\.e' '\.o' '\.pe' '\.po'      ! to remove pbs files \
                        "
convert.py_2to3_nb  =   "to convert ipynb files of python2 to python3\
                        \n\t\t to change .py file, use 2to3 in anaconda: source anaconda.sh and activate\
                        "

command.command     = "show Recent Command"
command.web_load    = "To load in web, copy files to ~/public_html/"

### removed
dirjob.dir_reset="reset dir as initial state by job: -j ai"
dirjob.dir_run="scan dir and run the same command for all the files such as\n\t\tqcout_mol_in.pl"

#classobj_dict={'EDIR':exe_dir, "MDIR":mod_dir}

def classify(Lclassify, work):
    
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
        if not Lclassify:
            for f in sort_exe:
                print(f"    {f}")
        else:
            
            for instance in MyClass.instances:
                #print(globals().keys())
                for gkey in globals().keys():
                    if gkey == instance.name:    # 'sge'(class instance) == 'sge'(string)
                        break
                if work != instance.name:
                    ckeys = dir_classify_n(sort_exe, instance.name, globals()[gkey], Lwrite=1)
                else:
                    ckeys = dir_classify_n(sort_exe, instance.name, globals()[gkey], Lwrite=0)
                    for ckey in ckeys:
                        print(f"    {ckey}.py[sh]\t:: {globals()[gkey].__dict__[ckey]}")
            print("  Remainder == ")
            for f in sort_exe:
                print(f"    {f}")
    if sort_mod:
        print("Module:: ")
        if not Lclassify:
            for f in sort_mod:
                print(f"    {f}")
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
            print("  Remainder == ")
            for f in sort_mod:
                print(f"    {f}")
    
    print("Instances:: ", end='')
    for instance in MyClass.instances:
        print(f"{instance.name}", end=' ')
    print("\n\t    -w for detail")
    print(f"#Comment: -c    for classification")
    '''
    if not spec == None:
        print(f"Detail for {spec}:: ")
        name_class = classobj_dict[spec]
        for key in name_class.__dict__.keys():
            print(f"    {key}\t:: {name_class.__dict__[key]}")
    '''
    return 0

def main():
    parser = argparse.ArgumentParser(description="display Usage for $SB/py_qcmo  ")
    parser.add_argument('-c', '--classify', action="store_true", help="classify files ")
    parser.add_argument('-w','--work',  help="several explanation option ")
    args = parser.parse_args()

    classify(args.classify,args.work)

if __name__ == "__main__":
    main()
