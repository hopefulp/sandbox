#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_files, MyClass

usage = {   'xyz22mol' : ' a.xyz[a.mol]\t# makes a.mol[a.xyz]',
            'qout_geo' : ' a.out\t# extract the last optimized geometry as a.xyz',
            'qout_non_opt': ' a.out\t# extract the last non-optimized geometry as a.xyz'
        }

qcout = MyClass()
qcout.qcout_geo="get optimized geometry (xyz, mol) from qc.out"
qcout.qcout_geo_nonopt="obtain last geometry when optimzation failed from qc.out"
qcout.qcout_mol_in="obtain input file from qcout file by option r=fname m=fname i=inputf"
xyz = MyClass()
xyz.xyz22mol="convert a.xyz to a.mol vice versa"
xyz.xyz2inp="convert a.xyz to a.inp"
xyz.xyz_angle="calculate angle"
xyz.xyz_dist="calculate distance"

classobj_dict={'XYZ':xyz, 'QC-OUT':qcout}


def if_usage(f):
    f_pre = f.split('.')[0]
    if f_pre in usage.keys():
        print("    {}".format(f), usage[f_pre])
    else:
        print("    {}".format(f))
def jobs(job,cclass,ifile,np):
    if job == None or re.search("cl", job):
        mdir = os.path.dirname(__file__)
        print(f"List directory of{mdir}")
        exe, mod = dir_files(mdir)
        print("Executable:: ")
        sort_exe = sorted(exe)
        sort_mod = sorted(mod)
        if job == None:
            for f in sort_exe:
                print("    {}".format(f))
        else:
            lxyz=[]
            lqcout=[]
            for f in sort_exe:
                if re.match('xyz',f):
                    lxyz.append(f)
                    continue
                elif re.match('qcout',f):
                    lqcout.append(f)
                    continue
                print("    {}".format(f))
            ### classify xyz files
            print("  {:<10}::".format('XYZ format'))
            for f in lxyz:
                if_usage(f)
            ### classify qout files
            print("  {:<10}::".format('QC-OUT'))
            for f in lqcout:
                if_usage(f)
        print("Module:: ")
        for f in sort_mod:
            print("    {}".format(f))
        if not cclass == None:
            print(f"Detail for {cclass}::")
            name_class = classobj_dict[cclass]
            for key in name_class.__dict__.keys():
                print(f"    {key}\t:: {name_class.__dict__[key]}")

        print("#Comment: -j run for 'how to run'\n\t  -j classify\n\t  -c class for detail ")
    elif job == 'run':
        com_serial="qchem {0}.in {0}.out &".format(ifile)
        com_parallel="mpirun -np {1} $QC/exe/qcprog {0}.in $QCSCRATCH > {0}.out &".format(ifile,np)
        print("{:^8}::".format('Chi'))
        print("\t{:^8}::".format('serial'), com_serial)
        print("\t{:^8}::".format('parallel'), com_parallel)

    

def main():

    parser = argparse.ArgumentParser(description="explanation for /pyqchem ")
    parser.add_argument('-j','--job',  help="qchem run in chi, mlet ")
    parser.add_argument('-c','--cname',  help="detail for each class ")
    parser.add_argument('-f','--infile',  help="qchem input file")
    parser.add_argument('-np','--nprocess', default=2, type=int, help="number of parallel process")
    args = parser.parse_args()

    jobs(args.job,args.cname,args.infile, args.nprocess)

if __name__ == "__main__":
    main()
