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
amp     = MyClass('amp')
qchem   = MyClass('qchem')
fconv   = MyClass('fconv')
general = MyClass('general')

amp.amp_run             ="amp_run.py -f OUTCAR -j tr -des gs -tef -nc 4 -hl 8 8 -el 0.001 0.003 -fl 0.00 -nt 4000 -ntr 1500 -dtype int -dl 3000 3500\
                        \n\t\t\t   amp_run.py -f OUTCAR -j te -tef -nc 1 -nt 4000 -ntr 100 -dtype int -dl 0 10\
                        \n\t\t\t    run amp for making db (gen descriptor), training, test, md etc\
                        \n\t\t\t    mem:\
                        \n\t\t\t\t2G for db\
                        \n\t\t\t\t12G for training: HL, total_images\
                        \n\t\t\t    -j tr -des gs -tef:\
                        \n\t\t\t\tuse train for gs test (in the script) and test force\
                        \n\t\t\t\ttr w. 100 images is enough?\
                        "
amp.amp_loop            ="amp_loop.py\
                        \n\t\t\t    loop for many situation used in SGE"
amp.amp_plot            ="amp_plot.py\
                        \n\t\t\t    plot amp_run.py test"
amp.amp_descriptor      ="called by amp_run.py\
                        \n\t\t\tprovide symmetry function to test diverse descriptor\
                        "
amp.amp_mod             ="called by amp_run.py\
                        \n\t\t\tprovide some functions\
                        "
general.make_dir        ="make_dir.py new_dir [old_dir] -w amp pbs -j tr\
                        \n\t\t\tto make new dir and copy or ln -s files\
                        \n\t\t\tUsage:\
                        \n\t\t\t    amp\
                        \n\t\t\t\tnew_dir -w amp -j ['tr','db','des','md']\
                        \n\t\t\t\t    des links OUTCAR\
                        \n\t\t\t\t    tr, db, md links OUTCAR, amp.db\
                        "

fconv.fconv2extxyz="convert file format to extxyz: im_format"
fconv.im2extxyz="convert IM file format to extxyz"
fconv.NucCarts2xyz="Not Used: use when View.xyz is not provide in Q-Chem, AIMD calculation"
fconv.aimd2extxyz="from Q-Chem, AIMD calculation with result of EComponents to extxyz format"

qchem.aimd = "AIMD run:: (chi::parallel) mpirun -np 4 $QC/exe/qcprog water1_aimd.in $QCSCRATCH/w1K400 > w1k400.out\
                    \n\t\t use high T to get proper Ek for small system, generate as the number of scratch folders as -np\
                    \n\t\t EComponent:: $1=time; $2=Etot; $3=Epot; $13=Ekin \
                    \n\t\t NR==1 {next} for skip\
                    \n\t\t NR==2 time=0.0 $3=Epot0(Epot=Etot of opt); $13=Ekin from T=Etot from Epot0\
                    \n\t\t awk '{ if(NR==1) {next} else if(NR==2) {epot0=$3; printf \"total energy %.7f\\n\", $13*2600} else {printf \"Epot %.7f Ekin %.7f Etot %.7f\\n\", ($3-epot0)*2600, $13*2600, ($2-epot0)*2600}}' EComponents\
                    "
qchem.NucCarts2xyz="Convert NucCarts (AIMD) to xyz format\
                    \n\t\t\tNucCarts2xyz.py -d dirname -a atom_series such as O H H"


classobj_dict={'AMP_RUN': amp, 'FILE_CONV': fconv, 'QCHEM': qchem} 

def fconvert_eg():
    print("    INF file to extxyz\n\t{} xxx.inf -a atom_list -y_bar [1,2,3]".format(amp_collection['file_conv']))
    print("\t: atom_list comes from 'molecules.py'")
    print("\t: -y_bar is options: 1 for energy only")
    return


def file_conversion():
    print("file conversion: {}".format(amp_collection['file_conv']))
    fconvert_eg()
    return 0


def run_amp(fname,HL,elimit,nc,Lgraph):
    hl = " ".join(HL)
    print("amp run        : {}".format(amp_collection['amp_run']))
    print("    For sample profile::\n\t{} {} profile".format(amp_collection['amp_run'],models['Diss_CHO']))
    print(f"    For job=tr(ain)::\n\t{amp_collection['amp_run']} {fname} tr -hl {hl} -el {elimit} -n 5 -nc {nc}")
    if Lgraph:
        print(f"    For job=te(st)::\n\t{amp_collection['amp_run']} {fname} te -hl {hl} -el {elimit} -n 5 -nc {nc} +g")
    else:
        print(f"    For job=te(st)::\n\t{amp_collection['amp_run']} {fname} te -hl {hl} -el {elimit} -n 5 -nc {nc} -g")
    print("    For job=md::\n\t{} {} md".format(amp_collection['amp_run'],models['ethylene']))
    print("    For validation::\n\t{} -h\n\t{} {} '4 4 4' 0.001 +g | sh".format(amp_collection['amp_valid'],amp_collection['amp_valid'],models['ethylene']))
    print("    For validation scan::")
    print("                    [scan|not] job-type[val|tr] fname n_core")
    print("\t{} -h\n\t{} scan val {} 6 ".format(amp_collection['amp_scan'],amp_collection['amp_scan'],models['water']))
    return 0

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
        for f in sort_dir:
            print(f"    {f}")
    
    if sort_exe: 
        print("Executable:: ")
        if not Lclassify:
            for f in sort_exe:
                print("    {}".format(f))
        else:
            ### confer "ini_pypbs.py", MyClass.instances is a list of string as class variables
            for instance in MyClass.instances:
                ### globals() includes MyClass() instances as keys
                for gkey in globals().keys():
                    if gkey == instance.name:
                        break
                if work != instance.name:
                    ckeys = dir_classify_n(sort_exe, instance.name, globals()[gkey], Lwrite=1) # globals()[instance.name] is not working
                else:
                    ckeys = dir_classify_n(sort_exe, instance.name, globals()[gkey], Lwrite=0)
                    for ckey in ckeys:
                        print(f"    {ckey}.py[sh]\t:: {globals()[gkey].__dict__[ckey]}")
            print("  == not analyzed")
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
            print("  == not analized ")
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
    print(f"#Comment: -c    for classification")

    return 0        

def main():

    parser = argparse.ArgumentParser(description="display Usage for ~/py_ai")
    parser.add_argument('-c', '--classify', action="store_true", help="classify files ")
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
