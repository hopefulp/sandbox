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
from common import list2dict, whereami


### INCAR ORDER for display
ordered_incar_keys=['SYSTEM','GGA','GGA_COMPACT','PREC','ALGO','NPAR','NCORE','NSIM','LPLANE','ISTART','ICHARG','ISPIN','ENCUT','NELM','NELMIN','NELMDL','EDIFF','ISYM','ADDGRID','LREAL','LASPH','LMAXMIX','NELECT','MAGMOM','NUPDOWN','ISMEAR','SIGMA','AMIX','BMIX','AMIN','IWAVPRE','ISIF','IBRION','NSW','POTIM','EDIFFG','TEBEG', 'TEEND','SYMPREC', 'SMASS', 'MDALGO', 'NBLOCK', 'NWRITE','LPETIM','LWAVE','LCHARG','LAECHG','LVTOT','LVHAR','LORBIT','NEDOS','EMIN','EMAX','LPARD','NBMOD','EINT','LSEPB','LSEPK','NFREE','LEPSILON','LMONO','IDIPOL','LDIPOL','GGA_COMPAT','LSORBIT','IVDW','LVDWSCS','LDAU','LDAUTYPE','LDAUL','LDAUU','LDAUJ','LDAUPRINT', 'ICORELEVEL', 'CLNT', 'CLN', 'CLL', 'CLZ', 'LSCALAPACK', 'IMAGES', 'SPRING', 'LCLIMB' ]

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
opt_active  = {'ISIF': 2, 'IBRION': 2, 'POTIM': 0.3, 'NSW': 1000, 'EDIFFG': -0.01}
### chg, sp
chg_out     = ['ISIF', 'IBRION', 'EDIFFG', 'POTIM']
chg_change  = {'LCHARG': '.T.'}
sp_out      = chg_out
sp_change   = {'LCHARG': '.F.'}
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
### KISTI: param in follows param out to replace
kisti_out = ['NPAR']
kisti_in = {'NPAR': ['NCORE', 20]}


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

def extract_kv_inline(line):
    '''
    Return key, value pair in a line
    return  key-value
            None
    '''
    if not '=' in line:
        return None, None
    else:
        lst = line.strip().split()
        key = lst[0]
    ### find key in a line whether active or not 
    if '#' in key:
        key = key.replace('#', '')
    if key.upper() in ordered_incar_keys:
        return key.upper(), lst[2]
    else:
        return None, None

### modifying INCAR by dict or extract value by key
def modify_incar_bykv(incar, inc_kv, icout=None, outf='INCAR.mod', mode='m'):
    '''
    incar       input file of INCAR
    inc_kv      list or dict for INCAR key-value
    icout       keys to be commented out
    mode        m for modify INCAR, inc_kv is dict
                e to extract value, inc_kv is list
    return      m lines - list of lines
                e list of values
    '''
    if mode == 'm':
    ### inc_kv should be dict or even number of list elements
        if isinstance(ickv, list):
            kws = list2dict(inc_kv)
            print(f"{kws} in {whereami()} at {__file__}")
        elif isinstance(inc_kv, dict):
            kws = inc_kv
    ### to extract, kws is keys
    else:
        kws = inc_kv
    #print(f"beginning: kws {kws}")
    iline=0
    with open(incar) as f:
        lines = f.readlines()
    ### save to lines or values in a list
    newlist=[]
    ### line analysis for INCAR
    for line in lines:
        iline += 1
        #print(f"{line_key}: {kws.keys()}")
        line_key, line_value = extract_kv_inline(line)
        if mode == 'm':
            if line_key:
                if line_key in kws.keys():
                    line = f" {line_key}   =  {kws[line_key]}   ! change in mod_incar.py\n"
                ### if activated and in the list of delete
                elif line_key in icout:
                    line = f" #{line.rstrip()}    ! change in modify_incar_bykv\n"
            newlist.append(line)
        else:
            if line_key in kws:
                print(f"found input key and values {line_key}, {line_value}")
                newlist.append(line_value)
            
    ### write new file
    if mode == 'm':
        with open(outf, 'w') as f:
            for line in newlist:
                #print(line)
                f.write(line)
                return 0
    else:
        print(f"returns list of values: {newlist} ")
        return newlist


def modify_incar_byjob(incar, job, dic=None, opt='ac', suff=None):
    #os.system(f"cp {INCAR} INCAR")
    #line_change_dict('INCAR', vasp_job.zpe)
    print(f"running in {whereami()}:{__file__}")
    if suff:
        outf = f"INCAR.{suff}"
    else:
        outf  = f"INCARnew.{job}"
    ### common dict for modification
    if f'{job}_change' in globals():
        paramch   = eval(f'{job}_change')
    ### in case comment out: list
    if f'{job}_out' in globals():
        paramout  = eval(f'{job}_out')
    ### in case uncomment: dict
    if f'{job}_active' in globals():
        paramin   = eval(f'{job}_active')
    if f'{job}_in' in globals():
        paramrep  = eval(f'{job}_in')
    #print(f"setting: paramch {paramch} add {dic}")
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
    i=0
    iline=0
    with open(incar) as f:
        lines = f.readlines()
    print(f"write to {outf}")
    ### open output file and write line by line of input INCAR
    #print(paramin.keys())
    with open(outf, 'w') as f:    
        for line in lines:
            iline += 1
            lst = line.strip().split()
            if len(lst) == 0:
                f.write(line)
                continue
            else:
                first_item = lst[0]
            ### [1] check paramin
            ### param uncomment: when active and change: first activate and change
            if 'paramin' in locals() and paramin :
                for key in paramin.keys():
                    tag_match = False
                    if key in first_item:
                        if first_item == f'#{key}':
                            newline = line[1:]
                            print(f"paramin:{newline} i {i} {iline}")
                            i += 1
                        ### if option == 'ac', in activation, replace at the same time
                        if opt and 'c' in opt:
                            line = replace_line(paramin, key, job) + "\n"
                        tag_match = True    # used in change
                        ### as for 1 key appearance, just apply once
                        #continue
                ### only apply once then remove key
                #if tag_match == True:
                #    del paramin[key]
            ### [2]                
            ### param change
            if 'paramch' in locals() and paramch:
                for key in paramch.keys():
                    tag_match = False
                    ### replace finds only the first letter is active
                    if re.match(key, first_item):
                        line = replace_line(paramch, key, job) + "\n"
                        tag_match == True
                ### when remove a key, all the keys are not applied
                #if tag_match == True:
                #    del paramch[key]
            ### [3]
            ### param comment out
            tag_out = False
            if  'paramout' in locals() and paramout :
                for param in paramout: # this is list
                    if param in line:
                        #print(f"param out:{param} i {i} {iline}")
                        i += 1
                        line = comment_out_line(first_item, job)
                        tag_out = True
                #if tag_match == True:
                #    paramout.remove(param)
                    
                        ### param rep comes together with param out
                        if param in paramrep:
                            addline = add_line( paramrep, param, job )
                            f.write(addline)
            #print(f"{iline} line: {line}")
            f.write(line)

    return outf


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

