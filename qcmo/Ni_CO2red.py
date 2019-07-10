#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
import sys
from common import dir_all, MyClass, yes_or_no

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

def jobs(job, job_level, model, nfile, Llink, link_file, f3_option,Lnbo,nao_chg,nao_chgsum, Lrun):
    One_model=class_dic[model]
    model_name="m"+str(model)
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
        mfiles=" ".join([m1.f3,m2.f3,m3.f3,m4.f3,m5.f3,m6.f3])
        atoms=" ".join(nao_chg)
        chgsum=" ".join(nao_chgsum)
        com=f"mplot_mo.py -f {mfiles} --{job} -ncg {atoms} -ncs {chgsum}"
        if job_level > 1:
            com +=  " | grep SUM | awk '{print $4}'"
        print("-js 2 for sum, 3 for plot")
        if yes_or_no(com+"  ::will you run? (for more option: use -js [int>1]"):
            os.system(com)
    return 0

def main():
    parser = argparse.ArgumentParser(description="display Usage for /Ni-CO2_red \n\tUsage:: sys.$0 -m {model index} -nf {int} [-l|-lf]  ")
    parser.add_argument('-j','--job', choices=['moc','nbo'],  help="run for MOC|NBO")
    parser.add_argument('-js','--job_level', default=1, type=int, help="run for MOC|NBO with more extension")
    parser.add_argument('-m', '--model', default='1', choices=['1','2','3','4','5','6'], help="as for 1-PP, 2-PPP, 3-PNP, 4-PNP, 5-NiFe, 6-SM")
    parser.add_argument('-n', '--nfile', default=5, type=int, help="list input qcout file ")
    parser.add_argument('-l', '--link', action='store_true', help="write link option")
    parser.add_argument('-lf', '--link_file', action='store_true', help="input link_id data file") 
    parser.add_argument('-bar', '--f3_option', action='store_true', help="use different A-fragment") 
    parser.add_argument('-nbo', action='store_true', help="Do NBO analysis")
    #parser.add_argument('-chg', '--chg_atoms', default=["C1","O2","O3","Ni4","N7"], nargs="*", help="Do NBO analysis")
    parser.add_argument('-chg', '--chg_atoms', default=["C","O","O","Ni","N","N","Fe"], nargs="*", help="Do NBO analysis")
    parser.add_argument('--chg_sum', default=["C","O","O"], nargs="*", help="atoms for charge sum such as 'C O O' or 'N N'")

    parser.add_argument('-r', '--run', action='store_true', help="input link_id data file") 
    args = parser.parse_args()

    if not args.job:
        print("input -j [moc|nbo]")
        sys.exit(1)

    jobs(args.job,args.job_level,args.model,args.nfile,args.link,args.link_file,args.f3_option,args.nbo,args.chg_atoms,args.chg_sum,args.run)

if __name__ == "__main__":
    main()
