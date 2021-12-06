### module for INCAR modification
'''
job zpe
    with CHGCAR: 0 1, without: 0 2
    NSW = 1 : not compatable with ICHRGE 11
        NSW=0, IBRION = -1
    K at gamma
    !NPAR   will use default
    IBRION = 5, 6: takes long time
job band structure
    pre-opt: CHGCAR
    ISTART = 1; ICHARG = 11; NSW=0; IBRION=-1; LORBIT=11;
    Do not change: PREC, ALGO
    CHGCAR error: input NGX,Y,Z 
    prepare KPOINTS
job vdw
    change only IVDW = 12
'''

import argparse
import fileinput
import re
import os
import sys

###### DICT for each job
### job_mod for the existing value
### job_comment for comment out
### job_uncomment to make it active

### Band structure: change params and comment out
band_change = {'ISTART': 1, 'ICHARG': 11, 'NSW': 0, 'IBRION': -1,'LCHARG': '.F.'}
band_out    = ['POTIM', 'ISIF', 'EDIFFG'] # for comment out
band_active = {'LORBIT': 11}
### DOS
dos_change = {'ISTART': 1, 'ICHARG': 11, 'NSW': 0, 'IBRION': -1,'ALGO':'Normal', 'EDIFF':1E-5, 'LCHARG': '.F.'}
dos_out    = ['POTIM', 'ISIF', 'EDIFFG'] # for comment out
dos_active = {'LORBIT': 11, 'NEDOS': 4001, 'EMIN':-8, 'EMAX':6}
### VDW
vdw_active  = {'IVDW': 12}
noD_out     = ['IVDW']
### OPT
opt_change  = {'NSW': 1000}
opt_active  = {'ISIF': 2, 'IBRION': 2, 'POTIM': 0.3}
### Cell OPT
copt_change = {'NSW': 1000}
copt_active = {'ISIF': 3, 'IBRION': 2, 'POTIM': 0.3}
### zpe
zpe_change  = {'ICHARG': 1,'POTIM': 0.015, 'IBRION':5, 'NSW':1}
zpe_active  = {'NFREE': 2}
zpe_out     = ['NPAR']
### MAGMOM
#mag_active = {'ISPIN': 2}
mag_change = {'ISPIN': 2}
### KISTI
kisti_out = ['NPAR']
kisti_active = {'NCORE': 20}


def replace_line(dic, key, job=None):
    newline = f" {key} = {dic[key]}         ! in {job}"
    return newline

def comment_out_line(line, job):
    newline = f"#{line} + in {job} "+ "\n"
    return newline
'''
def dict_update(param, dic):
    if param in globals():
        paramch   = eval(f'{job}_change')
    return paramch
'''
def modify_incar(incar, job, dic=None, opt='ac', suff=None):
    #os.system(f"cp {INCAR} INCAR")
    #line_change_dict('INCAR', vasp_job.zpe)
    if suff:
        outf = f"INCAR.{suff}"
    else:
        outf  = "INCARnew."+job
    ### common dict for modification
    if f'{job}_change' in globals():
        paramch   = eval(f'{job}_change')
    ### in case comment out: list
    if f'{job}_out' in globals():
        paramout  = eval(f'{job}_out')
    ### in case uncomment: dict
    if f'{job}_active' in globals():
        paramin   = eval(f'{job}_active')
    #print(f"setting: paramch {paramch}")
    if dic:
        print(f"is this True {dic}")
        ### append params
        if 'a' in opt:
            if 'paramin' in locals():
                paramin.update(dic)
                print("extend dict")
            else:
                paramin = dic
        ### change params
        elif 'c' in opt:
            if 'paramch' in locals():
                paramch.update(dic)
            else:
                paramch = dic
        ### out params
        elif 'o' in opt:
            if 'paramout' in locals():
                paramout.update(dic)
            else:
                paramout = dic
    #print(f"param active {paramin} param change {paramch}")
    # print(f"param comment out {paramout}")
    with open(incar) as f:
        lines = f.readlines()
    print(f"write to {outf}")        
    with open(outf, 'w') as f:    
        for line in lines:
            mline = line.strip()
            ### param uncomment: when active and change: first activate and change
            if 'paramin' in locals() and paramin :
                for key in paramin.keys():
                    tag_match = False
                    if key in mline:
                        #print(f"mline:{mline}")
                        if mline[0] == '#':
                            line = mline[1:] +'\n'
                            #print(f"paramin:{line}")
                        ### if option == 'ac', in activation, replace at the same time
                        if opt and 'c' in opt:
                            line = replace_line(paramin, key, job) + "\n"
                        tag_match = True
                ### only apply once then remove key
                #if tag_match == True:
                #    del paramin[key]
                            
            ### param change
            if 'paramch' in locals() and paramch:
                for key in paramch.keys():
                    tag_match == False
                    ### replace finds only the first letter is active
                    if re.match(key, mline):
                        line = replace_line(paramch, key, job) + "\n"
                        tag_match == True
                ### when remove a key, all the keys are not applied
                #if tag_match == True:
                #    del paramch[key]
            ### param comment out
            if  'paramout' in locals() and paramout :
                for param in paramout:
                    tag_match = False
                    if param in line:
                        line = comment_out_line(mline, job)
                        tag_match = True
                #if tag_match == True:
                #    paramout.remove(param)
                    
            #print(f"{line}", end='')
            f.write(line)

    return outf


def main():
    parser = argparse.ArgumentParser(description='test for INCAR change')
    parser.add_argument('inf', help='input incar file or directory')
    parser.add_argument('job', choices=["dos","band","pchg","chg","md","cont","ini","zpe","mol","wav",'vdw','noD','opt','copt','mag','kisti'], help='job for VASP')
    parser.add_argument('-suf', '--suffix', default='test', help='change the input filename')
    args = parser.parse_args()

    if os.path.isfile(args.inf):
        f = args.inf
    elif os.path.isdir(args.inf):
        f = args.inf + "/INCAR"

    modify_incar(f, args.job, suff=args.suffix)

if __name__ == "__main__":
    main()

