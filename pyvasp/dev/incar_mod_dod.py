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


### dict for dict_job selection
incar_change= {'band':band_change,                  'opt':opt_change}       # dict of dict
incar_out   = {'band':band_out                                      }   # dict of list
incar_active= {'band':band_active, 'vdw':vdw_active,'opt':opt_active}

def replace_line(dic, key, job=None):
    newline = f" {key} = {dic[key]}         ! in {job}"
    return newline

def comment_out_line(line, job):
    newline = f"#{line} + in {job} "+ "\n"
    return newline

def modify_incar(incar, job, comout=None):
    #os.system(f"cp {INCAR} INCAR")
    #line_change_dict('INCAR', vasp_job.zpe)
    ofname  = "INCARnew."+job
    ### common dict for modification
    if job in incar_change.keys():
        paramch   = incar_change[job]
    ### in case comment out: list
    if job in incar_out.keys():
        paramout  = incar_out[job]
    ### in case uncomment: dict
    if job in incar_active.keys():
        paramin   = incar_active[job]
    #print(f"setting:{paramin}")
    with open(incar) as f:
        lines = f.readlines()
    with open(ofname, 'w') as f:    
        for line in lines:
            mline = line.strip()
            ### param change
            if 'paramch' in locals() and paramch:
                for key in paramch.keys():
                    if re.match(key, mline):
                        line = replace_line(paramch, key, job) + "\n"
            ### param comment out
            if  'paramout' in locals() and paramout :
                for param in paramout.keys():
                    if param in line:
                        line = comment_out_line(mline, job)
            ### param uncomment
            if 'paramin' in locals() and paramin :
                for key in paramin.keys():
                    if key in mline:
                        #print(f"mline:{mline}")
                        if mline[0] == '#':
                            line = mline[1:] +'\n'
                            #print(f"paramin:{line}")
            #print(f"{line}", end='')
            f.write(line)

    return ofname

