#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_all, MyClass, dir_classify_n, whereami


amp_collection = {
    'file_conv':    'im2extxyz.py'  ,
    'amp_run':      'amp_ene.py'    ,
    'amp_valid':    'amp_validation.sh',
    'amp_scan':     'amp_loop.sh'}

models= {
    'ethylene': 'Ethylene.extxyz',
    'Diss_CHO': 'Diss_H2COH.extxyz',
    'water'   : 'water128.extxyz'}
### these are included as key in globals()
make    = MyClass('make')
run     = MyClass('run')
clean   = MyClass('clean')
modify  = MyClass('modify')
ase     = MyClass('ase')

make.vas_make_ini   ="===== Start VASP:: Make directory with input files\
                    \n\twill make or copy POSCAR, POTCAR, KPOINTS, INCAR\
                    \n\tPOSCAR: \
                    \n\t    structure - POSCAR_dirname\
                    \n\tPOTCAR: \
                    \n\t    later produced in the work directory by 'genpotcar.py -pp pbe'\
                    \n\tKPOINTS:\
                    \n\t    constructed by raw\
                    \n\tINCAR:\
                    \n\t    need to be prepared in advance\
                    \n\t    MAGMOM can be modified by: sed -i 's/.*MAGMOM.*/MAGMOM = 18*0 27*2 100*0/' INCAR\
                    \n\tcd dirname then make potcar (platinum) by run 'genpotcar.py -pp pbe'\
                    "
make.myvasp         ="modules for make vasp directory\
                    \n\tcalled by 'vas_make_ini.py'\
                    \n\tfunctions:\
                    \n\t    get_vasp_repository\
                    \n\t    make_kpoints\
                    \n\t    get_atoms_4pos(pos='POSCAR')\
                    \n\t\treturns atoms_list, Natoms_list\
                    \n\t    make_mag_4pos(poscar)\
                    \n\t\tcall get_atoms_4pos\
                    \n\t\treturn MAGMOM\
                    \n\t    make_incar\
                    \n\tRun by 'python -m myvasp -j getmag -p poscar'\
                    \n\t    to get MAGMOM\
                    "
make.vas_make_d2d   =" Make Vasp dir from the existing old dir\
                    \n\tvas_make_d2d.py old_dir new_dir job [options]\
                    \n\tjob = {ini, cont, md, ... post process}\
                    "

run.amp_env_run         ="amp_run.py in (envs) anaconda\
                        \n\t\t   when envs is not (base), detect envs and import proper module\
                        "
clean.clean         =" "
modify.pos_sort     ="pos_sort.py POSCAR\
                    \n\tsort atoms in POSCAR\
                    \n\treturns POSCARnew\
                    \n\t:when generate POSCAR via ASE w increasing supercell, atoms in order are replicated\
                    "
ase.ase_fconvert    =""
ase.ase_vasp        =""
ase.ase_zpe         =""
classobj_dict={'MAKE': make, 'RUN': run, 'CLEAN': clean} 

def classify(Lclassify, work, class_name, job, fname,HL, elimit, nc, Lgraph):
    
    mdir = os.path.dirname(__file__)
    print(f"List directory of {mdir} ")
    #exe, mod = dir_files(mdir)
    exe, mod, dirs, d_link = dir_all(mdir)
    sort_exe = sorted(exe)
    sort_mod = sorted(mod)
    sort_dir = sorted(dirs)

    if sort_dir:
        print("Directories:: ")
        if not Lclassify:
            for f in sort_dir:
                print(f"    {f}")
        else:
            for instance in MyClass.instances:
                for gkey in globals().keys():
                    if gkey == instance.name:
                        break
                if work != instance.name:
                    ckeys = dir_classify_n(sort_dir, instance.name, globals()[gkey], Lwrite=1) # globals()[instance.name] is not working
                else:
                    ckeys = dir_classify_n(sort_dir, instance.name, globals()[gkey], Lwrite=0)
                    for ckey in ckeys:
                        print(f"    {ckey}.py\t:: {globals()[gkey].__dict__[ckey]}")
            print("  == not classified")
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

    if job == 'amp':
        print("For AMP::")
        file_conversion()
        run_amp(fname,HL, elimit, nc, Lgraph)
    ### print dictionary here
    if job in classobj_dict.keys():
        name_class = classobj_dict[job]
        for key in name_class.__dict__.keys():
            print(f" {job} :: {name_class.__dict__[key]}")

    print("\nClass Instances:: ", end='')
    for instance in MyClass.instances:
        print(f"{instance.name}", end=' ')
    print("\n\t    -w for detail")
    #print(f"#Comment: -c    for not classification")

    return 0        

def main():

    parser = argparse.ArgumentParser(description="display Usage for ~/py_ai")
    parser.add_argument('-c', '--classify', action="store_false", help="classify files ")
    parser.add_argument('-w','--work',  help="several explanation option ")
    parser.add_argument('-j','--job',  help="[val,train,test] ")
    parser.add_argument('-cn', '--cname', help="detail for each class ")
    parser.add_argument('-f','--file',  help="input energy data file ")
    parser.add_argument('-hl','--hidden_layer',nargs='*', default=['4','4','4'], help="list of number of Hidden Layer")
    parser.add_argument('-el','--energy_limit',default=0.001, type=float,  help="energy_limit for training")
    parser.add_argument('-nc','--ncore',  help="number of parallel process")
    parser.add_argument('-g','--graph', action='store_true',  help="draw graph or not")
    args = parser.parse_args()

    classify(args.classify, args.work, args.cname, args.job,args.file, args.hidden_layer, args.energy_limit,args.ncore, args.graph )

if __name__ == "__main__":
    main()
