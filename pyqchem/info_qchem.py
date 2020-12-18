#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_all, MyClass, dir_classify_n, whereami
#from info_pbs import qsub -- happen collision in common.py
usage = {   'xyz22mol' : ' a.xyz[a.mol]\t# makes a.mol[a.xyz]',
            'qout_geo' : ' a.out\t# extract the last optimized geometry as a.xyz',
            'qout_non_opt': ' a.out\t# extract the last non-optimized geometry as a.xyz'
        }

qcout   = MyClass('qcout')
qcin    = MyClass('qcin')
xyz     = MyClass('xyz')
qcrun   = MyClass('qcrun')
install = MyClass('install')

qcout.qcout_geo="get optimized geometry (xyz, mol) from qc.out"
qcout.qcout_geo_nonopt="obtain last geometry when optimzation failed from qc.out"
qcout.qcout_nbo =   "Get NBO charges\
                    \n\tUsage::\
                    \n\t    qcout_nbo.py nbo -f qc.out -a atom_series -g atom_group_to_be_summed\
                    \n\t    qcout_nbo.py nbo -f PPN-Fe.out -a Fe P P N\
                    \n\tnbo\
                    \n\t-f one input file\
                    \n\t-a atom series\
                    \n\t-g atom group to be summed\
                    \n\timport analyze_nbo() in qcout_mod.py\
                    "
qcout.qcout_mod =   "Modules for qcout analysis\
                    \n\tNBO\
                    \n\t    analyze_nbo(file_pointer, atom_list, atom_group)\
                    \n\tMOC\
                    \n\t    imo_dic_basis_coeff(imoc_list)\
                    \n\t    trim_coeff(dic, n)\
                    \n\t    imo_basis(imoc_list, l_atoms, filter_tag, mo_id)\
                    \n\t    Cal_Ncore(dict)\
                    "
qcin.qcget_in   =   "obtain qchem input file from rem mol files\
                    \n\tOptions::\
                    \n\t    r=remfile\
                    \n\t    m=geometry file\
                    \n\t    i=output file, qchem input file of a.in\
                    \n\tUsage::\
                    \n\t    qcget_in.pl r=fname m=fname i=inputf\
                    "
qcin.qcget_georem = "soft link to 'qcget_in.pl'"                    
                
xyz.xyz22mol="convert a.xyz to a.mol vice versa"
xyz.xyz2inp="convert a.xyz to a.inp"
xyz.xyz_angle="calculate angle"
xyz.xyz_dist="calculate distance"

qcrun.qsub_server="rf. info_pbs.pbs.qsub_server.py\
                \n\tinfo_pbs.py -w qsub\
                "
qcrun.qcrun_mols    =   "Read *.mol files in directory and qsub\
                        /n/tUsage::\
                        /n/    qcrun_mols.sh (bash) [1:remfile]\
                        "
qcrun.qcrun_ins     =   "Read *.in in directory and qsub\
                        /n/tUsage::\
                        /n/    qcrun_ins.sh [1:
qcrun.aimd="refer to py_ai_ini.py"


install.gcc_se="\n\tPrerequisit:: \
                \n\t\tfftw in $QC_EXT_LIBS\
                \n\t\tOpenBLAS in $QC_EXT_LIBS for v.4 or /usr/local/lib64(root) for v.5\
                \n\t    Run  \"qchem -save a.in a.out savename\" works for only deprecate of $QCLOCALSCR"
install.gcc_pa="\n\tPrerequisit:: \
                \n\t\tfftw-mpi in $QC_EXT_LIBS\
                \n\t\tOpenBLAS-serial in /usr/local/lib64(root) for v.5\
                \n\t\topenmpi v.2.0.2\
                \n\t    Run  \"mpirun -np 4 a.in $QCSCRATCH > a.out\" remaines in $QCSCRATCH with -np directories"
install.intel_choi="locate properly\
                \n\t\tmight be slower than gcc, check it"
classobj_dict={'XYZ':xyz, 'QC-OUT':qcout}

def if_usage(f):
    f_pre = f.split('.')[0]
    if f_pre in usage.keys():
        print("    {}".format(f), usage[f_pre])
    else:
        print("    {}".format(f))
def classify(Lclassify, work, class_name, job, ifile,np):

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

    '''
    elif job == 'run':
        com_serial="qchem {0}.in {0}.out &".format(ifile)
        com_parallel="mpirun -np {1} $QC/exe/qcprog {0}.in $QCSCRATCH > {0}.out &".format(ifile,np)
        print("{:^8}::".format('Chi'))
        print("\t{:^8}::".format('serial'), com_serial)
        print("\t{:^8}::".format('parallel'), com_parallel)
    '''
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

    parser = argparse.ArgumentParser(description="explanation for /pyqchem ")
    parser.add_argument('-c', '--classify', action="store_false", help="classify files ")
    parser.add_argument('-w','--work',  help="several explanation option ")
    parser.add_argument('-j','--job',  help="qchem run in chi, mlet ")
    parser.add_argument('-f','--infile',  help="qchem input file")
    parser.add_argument('-np','--nprocess', default=2, type=int, help="number of parallel process")
    args = parser.parse_args()

    classify(args.classify,args.work, args.job,args.job,args.infile, args.nprocess)

if __name__ == "__main__":
    main()
