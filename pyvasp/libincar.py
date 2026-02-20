##!/home/joonho/anaconda3/bin/python
### module for INCAR modification
'''
    modify_incar_byjob
    modify_incar_bykv   new simple version to modify INCAR
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
import re
import os
import sys
from fline_edit import inf_print
from common import whereami
from libstr import li2dic
from copy import deepcopy

### INCAR ORDER for display
ordered_incar_keys=['SYSTEM','GGA','GGA_COMPACT','PREC','ALGO','NPAR','NCORE','NSIM','LPLANE','ISTART','ICHARG','ISPIN','ENCUT','NELM','NELMIN','NELMDL','EDIFF','ISYM','ADDGRID','LREAL','LASPH','LMAXMIX','NELECT','MAGMOM','NUPDOWN','ISMEAR','SIGMA','AMIX','BMIX','AMIN','IWAVPRE','ISIF','IBRION','NSW','POTIM','EDIFFG','TEBEG', 'TEEND','SYMPREC', 'SMASS', 'MDALGO', 'NBLOCK', 'NWRITE','LPETIM','LWAVE','LCHARG','LAECHG','LVTOT','LVHAR','LORBIT','NEDOS','EMIN','EMAX','LPARD','NBMOD','EINT','LSEPB','LSEPK','NFREE','LEPSILON','LMONO','IDIPOL','LDIPOL','GGA_COMPAT','LSORBIT','IVDW','LVDWSCS','LDAU','LDAUTYPE','LDAUL','LDAUU','LDAUJ','LDAUPRINT', 'ICORELEVEL', 'CLNT', 'CLN', 'CLL', 'CLZ', 'LSCALAPACK', 'IMAGES', 'SPRING', 'LCLIMB' ]

###### DICT for each job
### job_mod for the existing value
### job_comment for comment out
### job_uncomment to make it active

cont_change = {'ISTART': 1, 'ICHARG':0}     # read WAVECAR
### Band structure & DOS: change params and comment out
band_change = {'ISTART': 1, 'ICHARG': 11, 'NSW': 0, 'IBRION': -1,'LCHARG': '.F.'}
dosband_out    = ['POTIM', 'ISIF', 'EDIFFG'] # for comment out
band_active = {'LORBIT': 11}
dos_change = {'ISTART': 1, 'ICHARG': 11, 'NSW': 0, 'IBRION': -1,'ALGO':'Normal', 'EDIFF':1E-5, 'LCHARG': '.F.', 'ISMEAR': -5}
#dosband_out    = ['POTIM', 'ISIF', 'EDIFFG'] # for comment out
dos_active = {'LORBIT': 11, 'NEDOS': 1001, 'EMIN':-10, 'EMAX':10}
### VDW
vdw_active  = {'IVDW': 12}
noD_out     = ['IVDW']
### OPT
opt_change  = {'NSW': 1000}
opt_active  = {'ISIF': 2, 'IBRION': 2, 'POTIM': 0.3, 'NSW': 1000, 'EDIFFG': -0.01}
### CHG|SP: chg, pchg, pbchg [Bader] are different
### NSW = 0 for A+B, A, B
chg_out     = ['ISIF', 'IBRION', 'EDIFFG', 'POTIM']
chg_change  = {'LCHARG': '.T.'}
pchgB_out   =   ['ISTART']
pchgB_change = {'ICHARG': 11, 'LAECHG': '.TRUE.'}    
pchg_active = { 'LPARD': 'T', 'LSEPB' : 'FALSE', 'IBAND' : '304', 'LSEPK' : 'FALSE', 'KPUSE' : '1 2 3 4' } 
sp_out      = chg_out
sp_change   = {'LCHARG': '.F.'}
### spw to write PROCAR for matching kpoints in WAVCAR and PROCAR for pchg calculation
### for band calculation LWAVE = .TRUE. , LORBIT doesn't need
#spw_change = {'ISTART': 0, 'ICHARG': 2, 'LCHARG': '.TRUE.'}
spw_change = {'ISTART': 0, 'ICHARG': 2, 'LCHARG': '.TRUE.', 'LWAVE':'.TRUE.', 'NSW':0}

### Cell OPT
copt_change = {'NSW': 1000}
copt_active = {'ISIF': 3, 'IBRION': 2, 'POTIM': 0.3}
### zpe
zpe_change  = {'ICHARG': 1,'POTIM': 0.015, 'IBRION':5, 'NSW':1}
zpe_active  = {'NFREE': 2}
zpe_out     = ['NPAR']
### MAGMOM
#mag_active = {'ISPIN': 2}
mag_change = {'ISTART': 0, 'ICHARG': 1, 'ISPIN': 2, 'ISYM': 0}
### KISTI: param in follows param out to replace
kisti_out = ['NPAR']
kisti_in = {'NPAR': ['NCORE', 20]}
kisti_change = {'NCORE': 20}
### change: w.r.t. key -> value change
### active: if start with # -> remove #
### out   : comment out to remove the kv line
JOB_RULES = {
    "cont": {"change": cont_change,     "active": {},           "out": []},
    "dos" : {"change": dos_change ,     "active": dos_active,   "out": dosband_out},
    "band": {"change": band_change ,    "active": band_active,  "out": dosband_out},
    "mag":  {"change": mag_change ,     "active": {},           "out": []},
    "spw":  {"change": spw_change ,     "active": {},           "out": []},
    "kisti":{"change": kisti_change ,   "active": {},           "out": kisti_out},
}


def replace_line(dic, key, job=None):
    newline = f" {key} = {dic[key]}         ! in {job}"
    return newline

def add_line(dic, key, job=None):
    newline = f" {dic[key][0]} = {dic[key][1]}         ! in {job}\n"
    return newline
def comment_out_line(line, job):
    ''' if commented out already, return itself '''
    sline = line.strip()
    if re.match("#", sline):
        return line + "\n"
    else:
        newline = f"#{line} ! in {job} "+ "\n"
        return newline
'''
def dict_update(param, dic):
    if param in globals():
        paramch   = eval(f'{job}_change')
    return paramch
'''

comment = '!'

def extract_kv_inline(line):
    '''
    Return key, value, status
        if # commented, status returns 0, if not return 1
        if key not in ordered_incar_keys: return None
    '''

    linestrip=line.strip()
    status = 1  # key is active

    if not '=' in linestrip:
        return None, None, 0

    ### If '=' in line, cut '!'-after

    if comment in linestrip:
        keystring = linestrip.split(comment)[0]
    else:
        keystring = linestrip
    ### delete '#'
    if re.match('#', keystring):
        kvstring = keystring[1:].strip()
        status = 0
    else:
        kvstring = keystring
        status = 1
    #print(f"kvstring {kvstring}")

    ### Case '=', cut line before '!'
    ncount = kvstring.count('=')
    if ncount == 1:
        #if not ';' in line:
        lst = kvstring.split()   
        ### Case: key=value
        if '=' in lst[0]:
            kstr = re.split('=', lst[0])
            key = kstr[0]
            value = kstr[1]
        ### Case: key = value; 
        else:
            key = lst[0]
            value = lst[2]

        ### Case: only reserved keys are changed
        if key.upper() in ordered_incar_keys:
            return key.upper(), value, status
        else:
            print(f"{key} is not registered in ordered_incar_keys: line {line}")
            return None, None, 0
    else:
        if ncount == 2:
            print(f"Warning:: there are two k-v's in a line - {line}")
            pass
            #sys.exit(100)
            #print(f"two ='s in {linestrip}")
        return None, None, 0
        

### modifying INCAR by dict or extract value by key
def modify_incar_bykv(incar, inp_kv, icout=None, outf='INCAR.mod', mode='m'):
    '''
    incar       input file of INCAR
                INCAR need to have one key in a line
    inp_kv      list or dict for INCAR key-value
    icout       keys list to be commented out
    mode        m for modify INCAR, inp_kv is dict
                e,d to delete key, inp_kv is list
                if # key, remove #
    return      m output filename
                e,d list of values
    '''
    #inf_print(incar)
    if mode == 'm':
    ### inp_kv should be dict or even number of list elements
        if isinstance(inp_kv, list):
            kws = li2dic(inp_kv)
            print(f"{kws} in {whereami()} at {__file__}")
        elif isinstance(inp_kv, dict):
            kws = inp_kv
    ### to extract, kws is keys
    else:
        kws = inp_kv
    print(f"beginning: kws {kws}")
    iline=0
    with open(incar) as f:
        lines = f.readlines()
    ### save to lines or values in a list
    newlist=[]
    ### line analysis for INCAR
    all_keys=[]
    for line in lines:
        iline += 1
        ### if key is commented, activate
        #print(f"{line_key}: {kws.keys()}")
        line_key, line_value, status = extract_kv_inline(line)
        all_keys.append(line_key)
        if mode == 'm':
            if line_key:
                #print(f"line key {line_key}")
                if line_key in kws.keys():
                    print(f"line mod: {line_key} = {kws[line_key]}")
                    line = f" {line_key}   =  {kws[line_key]}   ! change in libincar.py\n"
                if icout:
                    if line_key in icout:
                        line = '#' + line
            ### Case: line has kws.keys(), it is modified
            newlist.append(line)
        else:
            if line_key in kws:
                #print(f"found input key and values {line_key}, {line_value}")
                newlist.append(line_value)
    print(f"{iline} was saved in newlist {len(newlist)}")
    ### write new file
    if mode == 'm':
        with open(outf, 'w') as f:
            for line in newlist:
                #print(line)
                f.write(line)
        print(f"{incar} was modified in {outf}")
        return outf
    else:
        print(f"returns list of values: {newlist} ")
        return newlist


def add_inckv_bysubjob(job, subjob, incdic):
    '''
    job md subjob: change INCAR.md by cool, heat, quenching
        quenching   NSW 1000
        cool, heat  default TEBIG, TEEND
    return INCAR dict to be modified by subjob
    '''
    if job == 'md':
        if subjob == 'quench':
            if not 'NSW' in incdic.keys():
                incdic['NSW'] = 1000
    
    return incdic

#def modify_incar_byjob(incar, job, outf='INCAR.new'):
def modify_incar_byjob(job):
    '''
    Modified    Do not write updated INCAR
                Return dict_change, incar_to_remove
    pass job with job + subjob
    job = cont, spw, spw2,
        spw     write WAVECAR, CHGCAR
    Return:
        dict_change (merged change+active)
        icout (list)
    '''

    if job not in JOB_RULES:
        print(f"no {job} in JOB_RULES in {whereami()}() in module {__file__}")
        sys.exit(101)

    rule = JOB_RULES[job]

    ### generate new dict (not touch original)
    dict_change = {}
    dict_change.update(rule.get("change", {}))
    dict_change.update(rule.get("active", {}))
    icout = list(rule.get("out", []))

    return dict_change, icout

def main():
    parser = argparse.ArgumentParser(description='test for INCAR change')
    parser.add_argument('inf', help='input incar file or directory')
    parser.add_argument('-j', '--job', choices=["dos","band","pchg","chg","md","cont","ini","zpe","mol","wav",'vdw','noD','opt','copt','mag','kisti'], help='job for VASP')
    parser.add_argument('-o', '--option', help='change the input filename')
    parser.add_argument('-kv', '--inc_dict', nargs='*', help='INCAR dict for modification')
    parser.add_argument('-suf', '--suffix', help='change the input filename')
    args = parser.parse_args()

    if os.path.isfile(args.inf):
        f = args.inf
    elif os.path.isdir(args.inf):
        f = args.inf + "/INCAR"

    if args.job:
        modify_incar_byjob(f, args.job, suff=args.suffix)
    else:
        modify_incar_bykv(f, args.inc_dict, mode=args.option)

if __name__ == "__main__":
    main()

