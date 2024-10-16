#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_files, MyClass, dir_classify


amp_collection = {
    'file_conv':    'im2extxyz.py'  ,
    'amp_run':      'amp_ene.py'    ,
    'amp_valid':    'amp_validation.sh',
    'amp_scan':     'amp_loop.sh'}

models= {
    'ethylene': 'Ethylene.extxyz',
    'Diss_CHO': 'Diss_H2COH.extxyz',
    'water'   : 'water128.extxyz'}

amp_run = MyClass()
amp_run.amp_ene="run amp for training, test, md etc"
amp_run.amp_loop="loop for many situation used in SGE"

fconv = MyClass()
fconv.fconv2extxyz="convert file format to extxyz: im_format"
fconv.im2extxyz="convert IM file format to extxyz"
fconv.NucCarts2xyz="Not Used: use when View.xyz is not provide in Q-Chem, AIMD calculation"
fconv.aimd2extxyz="from Q-Chem, AIMD calculation with result of EComponents to extxyz format"

qchem = MyClass()
qchem.aimd = "AIMD run:: (chi::parallel) mpirun -np 4 $QC/exe/qcprog water1_aimd.in $QCSCRATCH/w1K400 > w1k400.out\
                    \n\t\t use high T to get proper Ek for small system, generate as the number of scratch folders as -np\
                    \n\t\t EComponent:: $1=time; $2=Etot; $3=Epot; $13=Ekin \
                    \n\t\t NR==1 {next} for skip\
                    \n\t\t NR==2 time=0.0 $3=Epot0(Epot=Etot of opt); $13=Ekin from T=Etot from Epot0\
                    \n\t\t awk '{ if(NR==1) {next} else if(NR==2) {epot0=$3; printf \"total energy %.7f\\n\", $13*2600} else {printf \"Epot %.7f Ekin %.7f Etot %.7f\\n\", ($3-epot0)*2600, $13*2600, ($2-epot0)*2600}}' EComponents\
                    "
qchem.NucCarts2xyz="Convert NucCarts (AIMD) to xyz format\
                    \n\t\t\tNucCarts2xyz.py -d dirname -a atom_series such as O H H"


classobj_dict={'AMP_RUN': amp_run, 'FILE_CONV': fconv, 'QCHEM': qchem} 

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


def classify(Lclass, class_name, job, fname,HL, elimit, nc, Lgraph):
    
    mdir = os.path.dirname(__file__)
    print(f"List directory of {mdir} ")
    exe, mod = dir_files(mdir)
    sort_exe = sorted(exe)
    sort_mod = sorted(mod)
    print("Executable:: ")
    if not Lclass:
        for f in sort_exe:
            print("    {}".format(f))
    else:
        for key in classobj_dict.keys():
            fkey = dir_classify(sort_exe, key, classobj_dict, 1)
            if class_name:
                cname = classobj_dict[class_name]
                for key in cname.__dict__.keys():
                    if key in fkey:
                        print(f"    {key}.py\t:: {cname.__dict__[key]}")
        print("  == remainder ")
        for f in sort_exe:
            print(f"    {f}")
    if not sort_mod:
        print("Module:: not provided")
    else:
        print("Module:: ")
        if not Lclass:
            for f in sort_mod:
                print("    {}".format(f))
        else:
            for key in classobj_dict.keys():
                fkey = dir_classify(sort_mod, key, classobj_dict, 1)
                if class_name:
                    cname = classobj_dict[class_name]
                    for key in cname.__dict__.keys():
                        if key in fkey:
                            print(f"    {key}.py\t:: {cname.__dict__[key]}")
            print("  == remainder ")
            for f in sort_mod:
                print(f"    {f}")

    print("#Comment: try '-j amp'")
    if job == 'amp':
        print("For AMP::")
        file_conversion()
        run_amp(fname,HL, elimit, nc, Lgraph)
    ### print dictionary here
    if job in classobj_dict.keys():
        name_class = classobj_dict[job]
        for key in name_class.__dict__.keys():
            print(f" {job} :: {name_class.__dict__[key]}")

    if class_name == None:
        print(f"#Usage: -c for classification\
            \n\t-cn {classobj_dict.keys()} for detail ")
            

def main():

    parser = argparse.ArgumentParser(description="display Usage for ~/py_ai")
    parser.add_argument('-j','--job',  help="[val,train,test] ")
    parser.add_argument('-c', '--classify', action="store_true", help="classify files ")
    parser.add_argument('-cn', '--cname', help="detail for each class ")
    parser.add_argument('-f','--file',  help="input energy data file ")
    parser.add_argument('-hl','--hidden_layer',nargs='*', default=['4','4','4'], help="list of number of Hidden Layer")
    parser.add_argument('-el','--energy_limit',default=0.001, type=float,  help="energy_limit for training")
    parser.add_argument('-nc','--ncore',  help="number of parallel process")
    parser.add_argument('-g','--graph', action='store_true',  help="draw graph or not")
    args = parser.parse_args()

    classify(args.classify, args.cname, args.job,args.file, args.hidden_layer, args.energy_limit,args.ncore, args.graph )

if __name__ == "__main__":
    main()
