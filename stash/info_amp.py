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
ampplot = MyClass('ampplot')
ampga   = MyClass('ampga')
qchem   = MyClass('qchem')
fconv   = MyClass('fconv')
general = MyClass('general')
Gs      = MyClass('Gs')
aux     = MyClass('aux')
lammps  = MyClass('lammps')

ampga.gaamp               ="===== Genetic Algoritm for AMP =====\
                        \n\t\tclass GaAmp: instance runs for generation\
                        \n\t\t    amp_jobstr is required to run 'amp_run.py' with ga\
                        \n\t\tga.py: retains GA method\
                        \n\t\tqsub_restart.py for re-qsub\
                        "
ampga.ga_amp_run        ="\n\tga_amp_run.py -js sh -nc 10 -hl 5 -nn 5\
                        \n\t-js: job_submit: qsub, sh(node-bash)\
                        \n\t-nch: number of chromosome\
                        \n\t-hl: max number of hidden layers, -nn: max number of nodes in layer\
                        \n\t   old version of amp_run.py for GA\
                        "
amp.ampdir              ="run scripts by case\
                        \n\t\tUsage:: ampdir.sh $1\
                        \n\t\t$1::\
                        \n\t\t    wrapper: run amp_wrapper.py in dirs\
                        \n\t\t    wrapper_subdir: run amp_wrapper.py in subdirs\
                        "
amp.amp_wrapper         ="run 'amp_run.py' twice with job tr and te\
                        \n\t\t\timport 'amp_ini.py' for class Amp_string\
                        \n\t    Usage:: amp_wrapper.py -js [node|qsub] [-te] [-c] &\
                        \n\t\t-c: check amp_run.py string without job submit\
                        \n\t\t-te: run only 'te'\
                        \n\t\t: with waiting using time.sleep()\
                        n\t\t\t train: waiting file 'amp.amp' for sometimes the script is stopped when amp conv-failes\
                        "
amp.amp_ini             ="class Amp_string\
                        \n\t\t imported from 'amp_wrapper.py', 'ga_amp_run.py'\
                        "

amp.amp_run             ="amp running with 'training' or 'test' and 'make db', md, profile\
                        \n\tUsage::\
                        \n\t    (tr) amp_run.py -inf OUTCAR -j tr -nc 4 -hl 8 8 -el 0.001 -fl 0.1 0.04 SYM_FUNC DATA_SEL\
                        \n\t\tSYM_FUNC = '-des gs -pf log10 -pmod del -pmm 0.05 200.0 -pn 10'\
                        \n\t\t   '-des gs -pf powNN -pn 5'\
                        \n\t\tDATA_SEL = '-nt 4000 -ntr 1500 -dtype int -dl 1000 2500 3500 3600'\
                        \n\t\tmem:\
                        \n\t\t    2G for db\
                        \n\t\t    12G for training: HL, total_images\
                        \n\t    (te) amp_run.py -inf OUTCAR -j te DATA_SEL\
                        \n\t\tNot necessary: ncore, sym-function, \
                        \n\t    (md) amp_run.py -j md -inf w1.extxyz -i 800 -p amp.amp -dt 1 -ns 100\
                        \n\t\t: makes traj.traj file which can be read by $ase gui traj.traj\
                        \n\t\t  plot: 'showall.py -s -j amp -k md_anal'\
                        \n\tMODULE dependence :\
                        \n\t    after v9(latest before amp other module)\
                        \n\t    using venv in anaconda\
                        \n\t    (venv) python $SBamp/amp_run.py ... which overwrites shebang\
                        \n\t    e.g.: (ampG2off) python $SBamp/amp_run.py -inf OUTCAR -j tr -hl 4 -el 0.001 -fl 0.01 0.04 -nt 4000 -ntr 100 -dtype int\
                        "
amp.stat_check          ="Machine Learning: statistics calculation for script check\
                        \n\trun at work directory\
                        "
amp.amp_env_run         ="amp_run.py in (envs) anaconda\
                        \n\t\t   when envs is not (base), detect envs and import proper module\
                        "
amp.amp_loop            ="amp_loop.py\
                        \n\t\t\t\t: loop for many situation used in SGE"
amp.amp_descriptor      ="called by amp_run.py\
                        \n\t\t\tprovide symmetry function to test diverse descriptor\
                        "
amp.amp_anal            = "amp_anal.py -f ../OUTCAR -p hl44E0.001F0.1N100/amp-untrained-parameters.amp -im 1081 -ia 3 -t 'wrong F'\
                        \t\t   amp_anal.py -f ../OUTCAR -p hl44E0.001F0.1N100/amp-untrained-parameters.amp -im 1081 1083\
                        \n\t\t\t: plot fingerprint of an atom by: index of image, index of atom\
                        \n\t\t\t: without -ia atom index, fp ranges for all the kinds of atoms are plotted\
                        \n\t\t   amp_anal.py -im 0\
                        \n\t\t   amp_anal.py -im 0 1\
                        "
amp.amp_datamining      = "developed by 'amp_anal.py'\
                        \n\t\t\tamp_datamining.py -im 0  +g\
                        \n\t\t\tamp_datamining.py -im 0 1 +g\
                        \n\t\t\tamp_datamining.py -im 0 -ia 0 +g\
                        "

amp.amp_util             ="called by amp_run.py\
                        \n\t\t\tprovide some functions\
                        "
ampplot.ampplot_test    ="plot test file: 'test_energy.dat',  \
                        \n\truns 'myplot2D.draw_amp_twinx'\
                        \n\tUsage:\
                        \n\t    ampplot_test.py -f test_energy.dat -hl 10 10 -el 0.0001 -ntr 700 -nte 100 -tl '1 Water molecule'\
                        "
ampplot.ampplot_stat_dir=" plot amp test files: 'test_energy.pkl', 'test_force.pkl'\
                        \n\tUsage:\
                        "
ampplot.ampplot_dir     ="using input x-dir, y-[tr,te] for files\
                        \n\tOptions:\
                        \n\t    -p for directory prefix\
                        \n\t    -t title\
                        \n\t    -yd for multiple directories for multiplot\
                        \n\tUsage::\
                        \n\t    ampplot_dir.py -p NN -t 'Ndata 100 (Gs-pow)'\
                        \n\t    ampplot_dir.py -p NN -t 'Training Set: Gs-pow' -y tr -yd . Ndata300\
                        \n\t    ampplot_dir.py -p NN -t 'Test Set    : Gs-pow' -y te -yd . Ndata300\
                        "
Gs.amp_gsversion         =" === Gaussian Symmetry function for AMP descriptor ===\
                        \n\tamp_gversion.py to select different versions of gaussian.py in ~/amp\
                        "
Gs.my_descriptor        ="\n\tmodule for GS"

                        
general.make_dir        ="make_dir.py new_dir [old_dir] -w amp pbs -j tr\
                        \n\t\t\tto make new dir and copy or ln -s files\
                        \n\t\t\tUsage:\
                        \n\t\t\t    amp\
                        \n\t\t\t\tnew_dir -w amp -j ['tr','db','des','md']\
                        \n\t\t\t\t    des links OUTCAR\
                        \n\t\t\t\t    tr, db, md links OUTCAR, amp.db\
                        "
general.dir_scan        = "dir_scan.sh \
                        \n\t\t\tscan the present directory then, ampdb's\
                        \n\t\t\tcount the fingerprint files in .../loose\
                        \n\t\t\to check the calculated database\
                        "
lammps.generate_lammps  =   "    ====  LAMPHET    ====\
                            \n\tmake lammps input files: this is for water system only\
                            \n\t    INPUT:\
                            \n\t\tamp.amp - amp potential file\
                            \n\t\tstarting_configuration.traj - lammps trajectory file of one image\
                            \n\t    OUTPUT:\
                            \n\t\tpotential_Atomname - PROPhet pot file for lammps\
                            \n\t\tsystem.data - one of lammps input file: l.data, l.in\
                            \n\t    Linked to Lammps\
                            \n\t\tlmp_serial < system.in\
                            \n\t\t    : making lammps output file of log.lammps, dump.lammps\
                            \n\t\tthen use make_trajectory.py\
                            "
lammps.make_trajectory  =   "\n\t(modify script)Read lammps output and convert ASE trajectory\
                            \n\t    Usage:: make_trajectory.py\
                            "


fconv.fconv2extxyz="convert file format to extxyz: im_format"
fconv.im2extxyz="convert IM file format to extxyz"
fconv.NucCarts2xyz="Not Used: use when View.xyz is not provide in Q-Chem, AIMD calculation"
fconv.aimd2extxyz   =   "from Q-Chem, AIMD calculation with result of EComponents to extxyz format\
                        \n\t\t-f output_filename(default:test.extxyz) -d aimp_dir(can be parent dir)\
                        \n\t\ttry eV as it were without rescaling to minimum\
                        \n\t\tUsage::\
                        \n\t\t    aimd2extxyz.py -d water1_aimd\
                        \n\t\t    out: test.extxyz\
                        "
fconv.xyz2extxyz="from Q-Chem, AIMD change xyz to include force in extxyz format"



qchem.aimd = "AIMD run:: (chi::parallel) mpirun -np 4 $QC/exe/qcprog water1_aimd.in $QCSCRATCH/w1K400 > w1k400.out\
                    \n\t\t use high T to get proper Ek for small system, generate as the number of scratch folders as -np\
                    \n\t\t EComponent:: $1=time; $2=Etot; $3=Epot; $13=Ekin \
                    \n\t\t NR==1 {next} for skip\
                    \n\t\t NR==2 time=0.0 $3=Epot0(Epot=Etot of opt); $13=Ekin from T=Etot from Epot0\
                    \n\t\t awk '{ if(NR==1) {next} else if(NR==2) {epot0=$3; printf \"total energy %.7f\\n\", $13*2600} else {printf \"Epot %.7f Ekin %.7f Etot %.7f\\n\", ($3-epot0)*2600, $13*2600, ($2-epot0)*2600}}' EComponents\
                    "
qchem.NucCarts2xyz="Convert NucCarts (AIMD) to xyz format\
                    \n\t\t\tNucCarts2xyz.py -d dirname -a atom_series such as O H H"

aux.info_amp       ="info file of this directory"

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
