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
'''

import fileinput
import re
### band structure: change params and comment out

mod_band = {'ISTART': 1, 'ICHARG': 11, 'NSW': 0, 'IBRION': -1,'LCHARG': '.F.'}
commentout_band = ['POTIM', 'ISIF', 'EDIFFG'] # for comment out
add_band = {'LORBIT': 11}

incar_mod={'band':mod_band}       # dict of dict
incar_commentout={'band':commentout_band}   # dict of list
incar_add={'band':add_band}
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
    mod     = incar_mod[job]
    ### in case comment out: list
    if job in incar_commentout.keys():
        comout  = incar_commentout[job]
    ### in case uncomment: dict
    if job in incar_add.keys():
        add     = incar_add[job]

    with open(incar) as f:
        lines = f.readlines()
    with open(ofname, 'w') as f:    
        for line in lines:
            mline = line.strip()
            for key in mod.keys():
                if re.match(key, mline):
                    line = replace_line(mod, key, job) + "\n"
            if  'comout' in locals():
                for param in comout:
                    if param in line:
                        line = comment_out_line(mline, job)
            if 'add' in locals():
                for key in add.keys():
                    if key in line:
                        line = replace_line(add, key, job) + "\n"
            f.write(line)

    return ofname

