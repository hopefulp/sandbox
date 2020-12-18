#!/home/joonho/anaconda3/bin/python
''' modified fron Ni_CO2red.py'''
import argparse
import os
import re
import sys
from common import dir_all, MyClass, yes_or_no
from qcout_mod import *
import mystring

Ni_files = MyClass()

class Ni_CO2():
    atom_types_Ni=['Ni', 'Ni C1 O1 O2', 'C O1 O2']
    mo_types_Ni = ['ONE', 'SEL', 'ALL']
    atom_types_Fe=['Ni Fe', 'Ni Fe C1 O1 O2', 'C O1 O2']
    mo_types_Fe = ['SEL', 'SEL', 'ALL']
    def __init__(self, f1, f2, f3, f4, f5, A_type='Ni'):
        self.f1=f1
        self.f2=f2
        self.f3=f3
        self.f4=f4
        self.f5=f5
        self.nf1=f3
        self.nf3=" ".join([f2,f3,f4])
        self.nf3_=" ".join([f1,f3,f4])
        self.nf5=" ".join([f1,f2,f3,f4,f5])
        if A_type == 'Ni':
            self.atom_types=Ni_CO2.atom_types_Ni
            self.mo_types=Ni_CO2.mo_types_Ni
        elif A_type == 'Fe':
            self.atom_types=Ni_CO2.atom_types_Fe
            self.mo_types=Ni_CO2.mo_types_Fe

    def atom_type(self,nf):
        if nf == 1:
            self.atom_id="\""+self.atom_types[1]
            self.mo_type="\""+self.mo_types[1]
        elif nf == 3:
            self.atom_id="\""+"\" \"".join([self.atom_types[0],self.atom_types[1],self.atom_types[2]])+"\""
            self.mo_type="\""+"\" \"".join([self.mo_types[0],self.mo_types[1],self.mo_types[2]])+"\""
        elif nf == 5:
            self.atom_id="\""+ "\" \"".join([self.atom_types[0],self.atom_types[0],self.atom_types[1],self.atom_types[2],self.atom_types[2]])+"\""
            self.mo_type="\""+ "\" \"".join([self.mo_types[0],self.mo_types[0],self.mo_types[1],self.mo_types[2],self.mo_types[2]])+"\""
        return self.atom_id, self.mo_type 

m1=Ni_CO2("1-PP-A-relaxed.out", "1-PP-A.out", "1-PP.out", "1-PP-B.out", "CO2.out")
m2=Ni_CO2("2-PPP-A-relaxed.out", "2-PPP-A.out", "2-PPP.out", "2-PPP-B.out", "CO2.out")
m3=Ni_CO2("3-Me-PNP-A-relaxed.out", "3-Me-PNP-A.out", "3-Me-PNP.out", "3-Me-PNP-B.out", "CO2.out")
m4=Ni_CO2("4-PNP-A-relaxed.out", "4-PNP-A.out", "4-PNP.out", "4-PNP-B.out", "CO2.out")
m5=Ni_CO2("5-NiFe-A-sep-sp.out", "5-NiFe-A.out", "5-NiFe.out", "5-NiFe-B.out", "CO2.out", A_type='Fe')
m6=Ni_CO2("6-CC-NiFe-A-relax.out", "6-CC-NiFe-A.out", "6-CC-NiFe.out", "6-CC-NiFe-B.out", "CO2.out", A_type='Fe')

class_dic={"1":m1, "2":m2, "3":m3, "4":m4, "5":m5, "6":m6}


#models_dict

def nbos(job, job_level, files, nao_chg_atoms, nao_chgsum_atoms):
    if job == "moc":
        if nfile == 1:
            files = One_model.nf1

        elif nfile == 3:
            if not f3_option:
                files = One_model.nf3
            else:
                files = One_model.nf3_

        elif nfile == 5:
            files = One_model.nf5

        at_type, mo_type=One_model.atom_type(nfile)

        link = ""
        if Llink:
            link+=" -l "
        elif link_file:
            link+=f" -lf {model_name}-{nfile}f.dat"
        com=f"mplot_mo.py -f {files} -a {at_type} -t {mo_type} {link}"
        if yes_or_no(com+"  ::will you run?"):
            os.system(com)
    elif job == "nbo":
        mfiles=" ".join(files)
        atoms=" ".join(nao_chg_atoms)
        chgsum=" ".join(nao_chgsum_atoms)
        for f in files:
            with open(f, 'r') as inf:
                chgs = analyze_nbo(inf, nao_chg_atoms, nao_chgsum_atoms)
        print('\n'.join(chgs))                
        #com=f"mplot_mo.py -f {mfiles} --{job} -ncg {atoms} -ncs {chgsum}"
        if job_level > 1:
            com +=  " | grep SUM | awk '{print $4}'"
        print("-js 2 for sum, 3 for plot")
        #if yes_or_no(com+"  ::will you run? (for more option: use -js [int>1]"):
        #    os.system(com)
    return 0

def main():
    parser = argparse.ArgumentParser(description="Get NBO Charges ")
    parser.add_argument('job', default='nbo', choices=['moc','nbo'],  help="run for MOC|NBO")
    parser.add_argument('-js','--job_level', default=1, type=int, help="run for MOC|NBO with more extension")
    parser.add_argument('-f', '--files',  nargs='+', help="list input qcout file with nbo calculation ")
    #parser.add_argument('-chg', '--chg_atoms', default=["C1","O2","O3","Ni4","N7"], nargs="*", help="Do NBO analysis")
    parser.add_argument('-a', '--chg_atoms', default=["C","O","O","Ni","N","N","Fe"], nargs="*", help="Do NBO analysis")
    parser.add_argument('-g', '--chg_group', default=["C","O","O"], nargs="*", help="atoms for charge sum such as 'C O O' or 'N N'")

    args = parser.parse_args()

    if len(args.chg_atoms) == 1 and len(args.chg_atoms[0]) > 1:
        print("change atoms to atom list")
        chg_atoms = mystring.get_atomlist4str(args.chg_atoms[0])
    else:
        chg_atoms = args.chg_atoms

    print(chg_atoms)

    nbos(args.job,args.job_level,args.files,chg_atoms,args.chg_group)

if __name__ == "__main__":
    main()
