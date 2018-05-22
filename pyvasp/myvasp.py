### python-vasp module

import os
import json
import re

""" repository for many vasp script 
def get_vasp_repository():
def make_incar(dic, rw, iofile):
def make_kpoints(kp, MH|gamma, **args=[band])

"""

eps_H2O = 78.3


def get_vasp_repository():
    """ my vasp repository for POTCAR & KPOTINS """

    hostname = os.popen('hostname').read().rstrip()
    #print 'hostname is', hostname

    ### subprocess is not working
    #hostname = subprocess.popen('hostname', stdout=subprocess.PIPE, shell=True)
    #proc = subprocess.Popen(["cat", "/etc/services"], stdout=subprocess.PIPE, shell=True)
    #(out, err) = proc.communicate()
    #print "program output:", out
    if hostname == 'chi':
        ini_dvasp = '/Data/Bkaist/VaspINI'
    else:
        ini_dvasp = '/qcfs/joonho/VaspINI'

    print "vasp repository is ", ini_dvasp, ' in system ', hostname
    if not os.access(ini_dvasp, os.F_OK):
        print "Error:: the directory cannot be found\n stop"
        exit(1)
    return ini_dvasp

def make_kpoints(kp, method):
    """ Make KPOINTS file from """
    fname = 'KPOINTS'
    f = open(fname, 'w')

    if 'method' not in locals():
        method = "MH"

    f.write("Monkhorst Pack\n")                 # 1st line
    f.write("0\n")                              # 2nd line
    if re.match("M", method, re.IGNORECASE):
        f.write("Monkhosrt\n")                  # 3rd line
    else:
        f.write("Gamma\n")
    if len(kp) == 3:
        l = "  ".join(kp) + "\n"                # 4th liine
        f.write(l)  
    add = "000"                                 # 5th line
    l = "  ".join(add) + "\n"
    f.write(l)
    f.close()        
    return 0            

def make_incar(dic, rw, iofile):
    """ Make INCAR file from dictionary of command line arguments
        automatically write "incar.key"-json string
        when modify "incar.key (json string_", read iofile by --rw r
        filename can be change with -f, --iofile """

    ### args list:  1. nonmag, mag(fm), afm
    ###             2. ini, cont
    ###             3. precision
    ###             4. gga
    ###             5. dispersion
    ###             6.

    if rw == 'w':
        with open(iofile, 'w') as ofile:
            ofile.write(json.dumps(dic, indent=4))
    elif rw == 'r':
        dic = json.load(open(iofile))
    else:
        print "w/r error for iofile"
        exit(1)

    print dic

    fname = 'INCAR'
    f = open(fname, 'w')
    ### 0: system name
    comm = 'SYSTEM = ' + json.dumps(dic) + '\n\n'
    f.write(comm)
    ### 1: magnetism
    if dic['mag'] == 'nm':
        comm = '#MAGMOM:: natom*magnetism ... \n\n'
        f.write(comm)
    else:
        if dic['nmag']:
            natom = dic['nmag']
        else:
            question = "number of magnetic atom? "
            natom = int(str(raw_input(question)).strip())
        magmom = dic['magmom']
        comm = 'MAGMOM = '
        if dic['mag'] == 'fm':
            comm += str(natom) +'*'+magmom+' 999*0\n\n'
        else:
            for i in range(0, natom):
                if i % 2 == 0:
                    comm += str(magmom)+' '
                else:
                    comm += '-'+str(magmom)+' '
            comm += '999*0\n\n'                        
        f.write(comm)
    ### 2: start
    f.write('# continuation\n')
    if dic['ini'] == 'start':
        comm = 'ISTART = 0\nICHARG = 2\n'
    elif dic['ini'] == 'chg':
        comm = 'ISTART = 1\nICHARG = 1\n'
    else:       # cont goes to wav
        comm = 'ISTART = 1\nICHARG = 0\n'
            
    if dic['mag'] == 'nm':
        comm += 'ISPIN = 1\n\n'
    else:
        comm += 'ISPIN = 2\n\n'
    f.write(comm)

    ### 3: precision 
    f.write('# precision 1\n')
    com1 = 'ENCUT = %d\n' % dic['cutoff']
    com1 += 'PREC = %s\n' % dic['precision']
    com1 += 'ISMEAR = 0 ; SIGMA = 0.05\n'
    com1 += 'NELMIN = 4 ; NELM = 500\n'
    com1 += 'EDIFF = 1E-5 ; EDIFFG = -0.025\n\n'
    f.write(com1)
    f.write('# precision 2\n')
    com1 = 'GGA_COMPAT=.FALSE.\n'
    com1 += '#VOSKOWN = 1\n'
    com1 += 'ADDGRID = .TRUE.\n\n'
    f.write(com1)
    f.write('# mixing\n')
    com1 = 'LMAXMIX = 4\n'
    com1 += '#IMIX = 4; #AMIX = 0.2; #BMIX = 0.0001; #AMIX_MAG = 0.8; #BMIX_MAG = 0.0001\n\n'
    f.write(com1)
    f.write('# performance')
    com1 = 'ALGO = Fast\n'
    com1 += '#IALGO=48\n'
    com1 += 'NSIM = 4; NPAR = 4\n'
    com1 += 'LREAL = Auto; LPLANE = .TRUE.\n'
    com1 += 'LSCALAPACK = .FALSE.\n\n'
    f.write(com1)
    ###### 4: dft+D+U
    f.write('# functional (PE=PBE,RE=revPBE,b3=B3LYP,ML=vdw-df2, MK=rev-vdW-DF2, \n')
    f.write('# D correction if not vdW-DF (0-no, 1-d2, 11-d3_zero(Grimme), 12-de_BJ, 2-ts)\n')
    comm = 'GGA = '+dic['dft']+'\n'
    ### B3LYP
    if re.search('b3',dic['dft'],re.IGNORECASE):
        comm += 'LHFCALC = .TRUE.\n'
        comm += 'ALGO = D; TIME = 0.5\n'
        comm += 'AEXX = 0.2\n'
        comm += 'AGGAX = 0.72\n'
        comm += 'AGGAC = 0.81\n'
        comm += 'ALDAC = 0.19\n'
    elif re.search('hs',dic['dft'],re.IGNORECASE):
        comm += 'LHFCALC = .TRUE.\n'
        comm += 'ALGO = D; TIME = 0.4\n'
        comm += 'HFSCREEN = 0.207\n'
        comm += 'PRECFOCK = F\n'
        comm += 'NKRED = 2\n'
        comm += '#block the NPAR\n'
    elif re.search('e0',dic['dft'],re.IGNORECASE):
        comm += 'LHFCALC = .TRUE.\n'
        comm += 'ALGO = D; TIME = 0.4\n'
        comm += 'PRECFOCK = F\n'
        comm += 'NKRED = 2\n'
        comm += '#block the NPAR\n'
        
    if re.search('ML',dic['dft'],re.IGNORECASE) or re.search('MK',dic['dft'],re.IGNORECASE):
        comm += 'LUSE_VDW = .TRUE. \n'
        comm += 'Zab_vdW = 1.8867 \n'
        comm += 'AGGAC = 0.0000 \n'
        comm += 'LASPH = .TRUE. \n'
        if re.search('MK',dic['dft'],re.IGNORECASE):
            comm += 'PARAM1 = 0.1234 \n'
            comm += 'PARAM2 = 0.711357 \n'
    else:            
        if dic['dispersion']:
            comm += 'IVDW = 12    ! D2\n\n'
    comm += '\n'
    f.write(comm)
    ### DFT virtual orbital
    comm = '### DFT virtual orbital\n'
    comm += '#ALGO = Exact\n'
    comm += '#NBANDS = 64\n'
    comm += '#LOPTICS = .TRUE.\n\n'
    f.write(comm)

    f.write('### U-correction\n')
    if dic['uterm']:
        ldaj = 1.0
        ldau = ldaj + dic['uterm']
        comm = 'LDAU = .TRUE.\nLDAUTYPE = 2\nLDAUL = 2 -1 -1 -1\nLDAUU = '+str(ldau)+' 0 0 0\nLDAJ = 1 0 0 0\n\n'
    else:
        comm = '\n'
        f.write(comm)
    ### 5: Relaxation
    f.write('# Optimization\n')
    if dic['crelax'] == 'sp':
        comm = '#NSW = ; ISIF = ; IBRION = ; POTIM = \n\n'
    else:
        if dic['crelax'] == 'atom':
            isif = 2
        else:
            isif = 3
        comm = 'NSW = 999 ; ISIF = %d\n' % isif
        comm += 'IBRION = 2\nPOTIM = 0.3\n'
        comm += '#ADDGRID = .TRUE.\n\n'
    f.write(comm)
    ### Solvent effect
    f.write('# Solvent effect::higher Ecut is required (LSOL=.TRUE. EB_K for dielectric) \n')
    if dic['solvent']:
        comm = 'LSOL = .TRUE.\n'
        if dic['dielectric']:
            comm += 'EB_K = ' + dic['dielectric'] + '\n'
        comm += '#TAU = 0\n'
        comm += '#LRHOB = .TRUE\n\n'
        f.write(comm)
    ### AIMD
    f.write('# aimd - nvt:: nsw=5000; isif=0; ibrion=0;potim=2.0;tein=300;tebeg=300;teend=300;smass=0.05;isym=0\n\n')
    ###### 6: POST-SCF 
    f.write('###### POST SCF CALC :: SOC [DOS|PCHG|NEB]  #########\n')
    ### SOC
    f.write('### spin-orbit coupling   (LSORBIT=.TRUE.) \n')
    if dic['soc']:
        comm = 'LSORBIT = .TRUE.\n\n'
        f.write(comm)
    ### DOS
    if (dic['postscf'] == 'dos'):
        comm = '### DOS    ( EMIN, EMAX, NEDOS=6000) \n'
        comm += 'EMIN = -40\n'
        comm += 'EMAX = 20\n'
        comm += 'NEDOS = 6000\n'
    ### PCHG  
    elif (dic['postscf'] == 'pchg'):
        comm = '### Partial charge  (NBANDS, LPARD IBAND, KPUSE, LSEPB, LSEPK)\n'
        comm += 'NBANDS = 112\n'
        comm += 'LPARD = .TRUE.\n'
        comm += 'IBAND = 91 92 \n'
        comm += 'KPUSE = 1 2 3 4 \n'
        comm += 'LSEPB = .TRUE. \n'
        comm += 'LSEPK = .TRUE. \n'
    ### NEB
    elif (dic['postscf'] == 'neb'):
        f.write('### NEB \n')

    f.write(comm)

    ###### 7: LOGFILE 
    f.write('###### log\n')
    laechg = '.FALSE.'
    lwave = '.FALSE.'
    lcharg = '.FALSE.'
    if re.match('N',dic['log'],re.IGNORECASE):
        pass
    else:
        lwave = '.TRUE.'
        lcharg = '.TRUE.'
        if re.match('f', dic['log'], re.IGNORECASE):
            laechg = '.TRUE.'
    comm = 'LROBIT = 11\nLAECHG = %s\nLWAVE = %s\nLCHARG = %s\n\n' % (laechg, lwave, lcharg)
    f.write(comm)

    ###### Extras
    comm = '#IDIPOL = 3\n'
    f.write(comm)

    f.close()
    return 0
