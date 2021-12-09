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
check   = MyClass('check')

make.vas_make_ini   ="==================== Start VASP =======================================\
                    \n\t===== Make VASP initial directory ====\
                    \n\t===== Prepare POSCAR, POTCAR, KPOINTS, INCAR\
                    \n\tOPTIONS:\
                    \n\t    -s  poscar\
                    \n\t\tPOSCAR.dirname <- convention\
                    \n\t    -r  run qsub on any server\
                    \n\t    -o  qopt to use different qscript\
                    \n\t    -al stop questioning except -s poscar\
                    \n\tPOTCAR: \
                    \n\t    will be made by 'genpotcar.py -pp pbe' after cd and reading POSCAR\
                    \n\tKPOINTS:\
                    \n\t    constructed by raw\
                    \n\tINCAR:\
                    \n\t    need to be prepared in advance\
                    \n\te.g. (KISTI)\
                    \n\t    python $sbvas/vas_make_ini.py -r -o -al -s POSCAR.pdh2sidetop\
                    "
make.vas_make_d2d   =" Make Vasp dir from the existing old dir\
                    \n\tvas_make_d2d.py old_dir new_dir job [options]\
                    \n\tjob = {cont, md, ... post process}\
                    "
make.vas_make_incar ="\n\t    If not 'incar.key', make it, check it and modify it before run this again\
                    \n\t    based on 'incar.key', make INCAR\
                    "
make.vas_make_cont = " -d dir_list -j job -i incar_option -o option_poscar\
                    \n\tchange of d2d to make_cont\
                    \n\toptions:\
                    \n\t    -d all the directory list\
                    \n\t    -j [ini,cont] calls vasp_job_ini()\
                    \n\t    -j [incar,sp,chg,opt,copt,vdw,noD,mag,kisti] changes only INCAR vasp_job_incar()\
                    \n\t    -j [band,dos,zpe] changes INCAR, KPOINTS,POSCAR in vasp_jobs()\
                    \n\t\tonly INCAR modified: vdw,noD, opt,copt, mag, kisti,incar\
                    \n\t\t    incar: no preexist values, given from cli\
                    \n\t    -il change of INCAR list for delete\
                    \n\t    -id change of INCAR dict\
                    \n\tUsage:\
                    \n\t    pypath.sh vas_make_cont.py -d SnO2sc22FH -j band -i i\
                    \n\t    python $sbvas/vas_make_cont.py -d sc34 -nd sc34E5 -j incar -id '{\"ENCUT\": \"500\"}' -io c\
                    \n\tJobs Detail:\
                    \n\t    ini, cont: (def ini)\
                    \n\t\tcopy odir [INCAR, KPOINTS, POTCAR] w. given -s POSCAR\
                    \n\t\t     if not -s POSCAR, use odir/POSCAR(ini) or CONTCAR(cont)\
                    \n\t\t<eg> -d odir -j ini -s POSCAR.newmodel\
                    \n\t\t    : default ndir is newmodel\
                    \n\t    opt\
                    \n\t\tcopy odir [CONTCAR, KPOINTS, POTCAR] w. change of INCAR\
                    \n\t    ZPE:\
                    \n\t\tafter opt, calculage zpe\
                    \n\t    band:\
                    \n\t\tafter LCHARG=.TRUE., calculate band structure\
                    \n\tChange INCAR:\
                    \n\t    import mod_incar\
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
make.mod_incar     ="module for INCAR modification\
                    \n\tdef modify_incar(incar_in, job, dic=None, opt='ac')\
                    \n\t    incar_in is modified by job_dict\
                    "
make.mod_poscar    ="module for POSCAR modification\
                    \n\tget_atoms_poscar\
                    \n\tget_poscar(poscar)\
                    \n\t    : read input POSCAR.job and write to wdir as 'POSCAR'\
                    \n\tpos2dirname(poscar)\
                    \n\t    : input POSCAR.job gives dirname of 'job'\
                    \n\tfixedMD_POSCAR(poscar, atom, atoms=None)\
                    \n\t    : poscar is modified\
                    \n\t    : all the atoms except atom will be fixed for ZPE\
                    "
run.amp_env_run     ="amp_run.py in (envs) anaconda\
                    \n\t\t   when envs is not (base), detect envs and import proper module\
                    "
run.vas_qsub        = " run vasp in queue\
                    \n\t called by vas_make_cont.py\
                    "
run.envvasp         = " imported from vasp run, make\
                    "
clean.clean         =" "
modify.pos_sort     ="pos_sort.py POSCAR\
                    \n\tsort atoms in POSCAR\
                    \n\treturns POSCARnew\
                    \n\t:when generate POSCAR via ASE w increasing supercell, atoms in order are replicated\
                    "
check.incar_diff   ="diff_incar.py INCAR1 [INCAR2] -k keys -a -s\
                    \n\tOptions:\
                    \n\t    one dir : show INCAR\
                    \n\t    two dirs: compare two INCAR\
                    \n\t    -k  keys: to check whether the keys exist\
                    \n\t    -a  : search all directories as for -k keys\
                    \n\tUsage:\
                    \n\t    incar_diff.py FPtb2H2 FPtb2H2hb\
                    \n\t    incar_diff.py -a -k ENCUT ISTART\
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
