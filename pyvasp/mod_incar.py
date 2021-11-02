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

import fileinput
import re
import sys

###### DICT for each job
### job_mod for the existing value
### job_comment for comment out
### job_uncomment to make it active
### dict for band structure: change params and comment out
band_change = {'ISTART': 1, 'ICHARG': 11, 'NSW': 0, 'IBRION': -1,'LCHARG': '.F.'}
band_out    = ['POTIM', 'ISIF', 'EDIFFG'] # for comment out
band_active = {'LORBIT': 11}
### dict for vdw
vdw_active  = {'IVDW': 12}
### dict for opt
opt_change  = {'NSW': 1000}
opt_active  = {'ISIF': 2, 'IBRION': 2, 'POTIM': 0.3}
### dict for mag
#mag_active = {'ISPIN': 2}
mag_change = {'ISPIN': 2}


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
def modify_incar(incar, job, dic=None, opt='ac'):
    #os.system(f"cp {INCAR} INCAR")
    #line_change_dict('INCAR', vasp_job.zpe)
    ofname  = "INCARnew."+job
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
        if 'a' in opt:
            if 'paramin' in locals():
                paramin.update(dic)
                print("extend dict")
            else:
                paramin = dic
        elif 'c' in opt:
            if 'paramch' in locals():
                paramch.update(dic)
            else:
                paramch = dic
        elif 'o' in opt:
            if 'paramout' in locals():
                paramout.update(dic)
            else:
                paramout = dic
    #print(f"param active {paramin} param change {paramch}")
    # print(f"param comment out {paramout}")
    with open(incar) as f:
        lines = f.readlines()
    with open(ofname, 'w') as f:    
        for line in lines:
            mline = line.strip()
            ### param uncomment: when active and change: first activate and change
            if 'paramin' in locals() and paramin :
                for key in paramin.keys():
                    if key in mline:
                        #print(f"mline:{mline}")
                        if mline[0] == '#':
                            line = mline[1:] +'\n'
                            #print(f"paramin:{line}")
                        ### if option == 'ac', in activation, replace at the same time
                        if opt and 'c' in opt:
                            line = replace_line(paramin, key, job) + "\n"
                            
            ### param change
            if 'paramch' in locals() and paramch:
                for key in paramch.keys():
                    ### replace finds only the first letter is active
                    if re.match(key, mline):
                        line = replace_line(paramch, key, job) + "\n"
            ### param comment out
            if  'paramout' in locals() and paramout :
                for param in paramout:
                    if param in line:
                        line = comment_out_line(mline, job)
            #print(f"{line}", end='')
            f.write(line)

    return ofname

