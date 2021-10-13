#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_all, MyClass, dir_classify_n, whereami

comment = MyClass('comment')
dirjob  = MyClass('dirjob')
jobfile = MyClass('jobfile')
convert = MyClass('convert')
command = MyClass('command')
string  = MyClass('string')

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



dirjob.dir_clean      =   "=== DIRECTORY JOB ==\
                        \n\t\tclean one directory\
                        \n\t\t    by -prefix -suffix -middle match -e excluded -ef 'exclude these files' -work {qchem,ai,amp,pbs} -j rm[mv] -jd new_dir\
                        \n\t\tOptions::\
                        \n\t\t    -d input directory: default=pwd\
                        \n\t\t    -y execute command without asking: default-asking\
                        \n\t\t    -w multiple jobs in [amp|pbs|slurm|vasp|lammps] add more in case extension\
                        \n\t\t    -sw subwork: amp-ini,ag\
                        \n\t\t    -j [rm,mv,cp,ln] default='rm'\
                        \n\t\t    -a: remove all such as 'rm -r'\
                        \n\t\t    -ef: excluded files\
                        \n\t\t\tln: in case the change of dirname, link files are broken\
                        \n\t\tUsage::\
                        \n\t\t    dir_clean.py [dir] -w amp -sw ini\
                        \n\t\t    dir_clean.py -d NN20 -w amp -j ln -y\
                        \n\t\t    dir_clean.py -w vasp -ef CHGCAR\
                        \n\t\t    (?) clean1d.py -s out -ef 6-CC-NiFe-A-relax.out 5-FePNP-CO2.out -j mv -jd j631gs_v3.2\
                        "
dirjob.clean_dirs   =   "clean dirs:: same with 'clean1d.py'\
                        \n\t\t input directory is list\
                        "
dirjob.clean_recur  = "clean one dir recursively\
                        \n\t\timport dir_clean() from clean1d.py to clean one directory\
                        \n\t\tUsage::\
                        \n\t\t    clean_recur.py -w pwd\
                        \n\t\t    clean_recur.py -d NN20 -w amp -j ln -y\
                        \n\t\t\trf. clean1d.py for options\
                        "
dirjob.cli_dir      =   "simple command inside directory\
                        \n\t\tdir1_cli.sh gitpush\
                        "
dirjob.git_push     =   "cli_dir.py gitpush was copied for only gitpush\
                        \n\t\t run in git directory\
                        "
dirjob.cli_dirs      =   "basic command line interface for all files in directory\
                        \n\t\tmodify script for all the files/selected files\
                        \n\t\tdir_cli.sh [0:vmake 1:incar 2:qsub 3:cp 4:rm 5:chmod\
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

string.common       =   "module for dir, string\
                        \n\t    Classes\
                        \n\t\tMyClass_obj: class has its name as class attribute\
                        \n\t\tMyClass(MyClass_obj): inherits MyClass_obj\
                        \n\t\tMyClass_str(dict): is not working\
                        \n\t    Functions::\
                        \n\t\tsearch_dirs(dir_prefix, filename)\
                        \n\t\tyes_or_no(string): get y/n from stdio\
                        \n\t      Get files from directory::\
                        \n\t\tget_files_type(filetype, dirname)\
                        \n\t\tget_files_prefix(prefix, dirname, Lshow, Ldir)\
                        \n\t      Filename::\
                        \n\t\tf_ext(fname): returns extension using [-1]\
                        \n\t\tf_root(fname): returns filename without extension\
                        \n\t\tfname_decom(fname): returns (fname, extension)\
                        \n\t\tetc\
                        "

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
    parser.add_argument('-c', '--classify', action="store_false", help="classify files ")
    parser.add_argument('-w','--work',  help="several explanation option ")
    args = parser.parse_args()

    classify(args.classify,args.work)

if __name__ == "__main__":
    main()
