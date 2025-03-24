#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_all, MyClass, dir_classify_n, whereami

comment = MyClass('comment')
dirjob  = MyClass('dirjob')
filejob = MyClass('filejob')
linejob = MyClass('linejob')
convert = MyClass('convert')
command = MyClass('command')
server  = MyClass('server')

server.server_env   =   "=== Server-related ===\
                        \n\t\tobtain server hostname, home etc\
                        "
server.pypath       =   "pypath.sh $python_command\
                        \n\t\tUsage::\
                        \n\t\t    $pypath.sh $(which script.py) args\
                        \n\t\tError::$script.py args -> bad interpreter: No such file due to shebang\
                        "
server.git_push     =   "cli_dir.py gitpush was copied for only gitpush\
                        \n\t\t run in git directory\
                        "

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
comment.info_common=   "\tInformation file for $SB/pycommon"


dirjob.dir_clean      =   "=============================== DIRECTORY JOB =====================================\
                        \n\t\tclean one directory by one option: works, prefixes, suffixes, matches\
                        \n\t\t    -work {qchem,ai,amp,pbs} -j rm[mv] -jd new_dir ? by -prefix -suffix -middle match -e excluded -ef 'exclude these files' \
                        \n\t\tOptions::\
                        \n\t\t    exclusive:\
                        \n\t\t\t-w multiple works in [amp|pbs|slurm|vasp|nc|lammps] add more in case extension\
                        \n\t\t\t-p multiple prefixes\
                        \n\t\t\t-s multiple sufixes\
                        \n\t\t\t-m multiple matches\
                        \n\t\t    -j [rm,mv,cp,ln] default='rm'\
                        \n\t\t    -id include dir such as 'rm -r'\
                        \n\t\t    -nd input directory: default=pwd\
                        \n\t\t    -y execute command without asking: default-asking\
                        \n\t\t    -sw subwork: amp-ini,ag\
                        \n\t\t    -a: remove all such as 'rm -r'\
                        \n\t\t    -ef: excluded files\
                        \n\t\t\tln: in case the change of dirname, link files are broken\
                        \n\t\tUsage::\
                        \n\t\t    dir_clean.py [dir] -w amp -sw ini\
                        \n\t\t    dir_clean.py -d NN20 -w amp -j ln -y\
                        \n\t\t    dir_clean.py -w vasp -ef CHGCAR\
                        \n\t\t    dir_clean.py -p d2709 -id\
                        \n\t\t    (?) clean1d.py -s out -ef 6-CC-NiFe-A-relax.out 5-FePNP-CO2.out -j mv -jd j631gs_v3.2\
                        \n\t\t    dir_clean.py -w vasp -d\
                        \n\t\t    : remove only vasp output files\
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
dirjob.dir_fname    =   "Treat Dir without work-style\
                        \n\t\tUsage:: dir_fname.py {ls, mv, rm, rename, cp, chmod} -options\
                        \n\t\t    mv    : mv to dir\
                        \n\t\t    rename: change filename\
                        \n\t\tOptions::\
                        \n\t\t    -[p|s|m] for matching type\
                        \n\t\t    -id to include dir in scanning dir\
                        \n\t\t    -v  inverse the matching\
                        \n\t\t    -r  recursive for subdirectroies\
                        \n\t\t    -st style=[ap:append, rp:replace, mo: mode\
                        \n\t\t    -rw replacement word, if None, replace becomes remove\
                        \n\t\t    -d, -nd  dirname for mv\
                        \n\t\t    -ip include_parents directory for matching is 'suffix'\
                        \n\t\t    -e exception list using matching\
                        \n\t\t    -eo default=m exception by matching or fullname\
                        \n\t\tUsage::\
                        \n\t\t    dir_fname.py rename -p G4 -st rp -rw G2 -id         ! rename with full replacement\
                        \n\t\t    dir_fname.py rename -p ToBeDelete_ -st rp -r        ! without -rw, -p is deleted in fname\
                        \n\t\t    dir_fname.py rename -p G4 -st ap -rw vdw -id        ! append new word after dir and fname\
                        \n\t\t    dir_fname.py rename -m sc34c -rp sc34 -id -e sc34ch ! rename dir & file with exception\
                        \n\t\t    dir_fname.py rm -m '\.e' '\.o' '\.pe' '\.po'        ! to remove pbs files \
                        \n\t\t    dir_fname.py mv -s .cont -id -d tmpcont -ip         ! to move include rootname dir\
                        \n\t\t    dir_fname.py mv -p HfSe2sc34 -nd HFSe2sc34 -id\
                        "
filejob.fline_edit  =   "Job to treat file:\
                        \n\tfind a line and substitute\
                        "
filejob.fline_sub   =   "template for line substitution"
filejob.fline_cut   =   "extract a certain part in a file\
                        \n\tusing keyword for start and end\
                        "
filejob.fline_part =    "extract a part from files: type=molden|band\
                        \n\tfname -j jobtype[molden, band] -k1 keyword -k2 keyword -i band_index\
                        \n\t    jobtype can be given by fname\
                        \n\t    keyword can be given by fname\
                        \n\tUsage::\
                        \n\t    fline_part.py BAND.dat -i 24\
                        \n\t\tto extract certain band from BAND.dat into 'BAND.idb01' only one band\
                        "                        
filejob.fline_shebang   =   "Change shebang to make .py executable"
filejob._extract_line   = "to extract a certain part in a file"
filejob.f_kw        =   "gather keywords to cut part of file\
                        \n\tQ-Chem outfile, also refer to ~/dev/\
                        \n\tBAND.dat from vaspkit to get part of bands\
                        "
filejob.fmath           = "to selectively modify column values\
                        \n\tUsage::\
                        \n\t    $fmath.py fname.dat -m prod -v 10 -xgt 0.7 \
                        \n\t\tmakes fname_new\{value\}.dat\
                        \n\tOptions::\
                        \n\t    -m:  [prod,add,div,sub]\
                        \n\t    -v:  value to be operated, default=5\
                        \n\t    -x: xcol to find region: default=0\
                        \n\t    -y: ycol to be modified: default=1\
                        \n\t    (Exclusive group)\
                        \n\t    -xij: xi and xj (nargs=2)  with values(float) or index(int)\
                        \n\t    -hlt: region x < 0.5\
                        \n\t    -hgt: regipn 0.5 < x\
                        \n\t    -xlt: value, region x < value\
                        \n\t    -xgt: value, region 0.5 < x\
                        \n\te.g.:\
                        \n\t    $fmath.py 2ldosa25_N72_192V-0.92.dat -m prod -v 10 -xgt 0.7 \
                        \n\t\tas for x values in col[0], if 0.7 < x, col[1]*3\
                        "

linejob.parsing     =   "line(string) parsing::\
                        \n\tdef is_int_el: \
                        \n\tdef is_int:\
                        \n\tdef is_there_char:\
                        \n\tdef str_decom:\
                        \n\tdef convert_s2l:\
                        \n\tdef get_atomlist4str:\
                        "
convert.py_2to3_nb  =   "to convert ipynb files of python2 to python3\
                        \n\t\t to change .py file, use 2to3 in anaconda: source anaconda.sh and activate\
                        "
convert.perl2python  =   " perl script to convert perl to python"
command.command     = "show Recent Command"
command.web_load    = "To load in web, copy files to ~/public_html/"
command.repeat_command = "to repeat cli command"

dirjob.common       =   "module for directory & string\
                        \n\t    Classes\
                        \n\t\tMyClass_obj: class has its name as class attribute\
                        \n\t\tMyClass(MyClass_obj): inherits MyClass_obj\
                        \n\t\tMyClass_str(dict): is not working\
                        \n\t    Functions::\
                        \n\t\tsearch_dirs(dir_prefix, filename)\
                        \n\t\tyes_or_no(string): get y/n from stdio\
                        \n\t    Get files from directory::\
                        \n\t\tget_files_type(filetype, dirname)\
                        \n\t\tget_files_prefix(prefix, dirname, Lshow, Ldir)\
                        \n\t    Filename::\
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
