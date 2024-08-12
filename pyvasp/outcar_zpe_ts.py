#!/home/joonho/anaconda3/bin/python

# Updated by Min Jong Noh : 2019.03.08 mailto:starnmj@kaist.ac.kr
# Updated by Joonho Park  : 2021.11.23, Fvib was modified, 

import os, sys, math
import argparse

#### Get a total number of lines in read-file ####

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i+1

def cal_zpe_ts(indir, natom):
    cwd = os.getcwd()
    if os.path.isdir(indir):
        os.chdir(indir)
        fname = "OUTCAR"
    elif os.path.isfile(indir):
        fname = indir

    #### THz Data mining in OUTCAR ####
    os.system(f'grep THz {fname} > THz.txt')

    TN_lines = file_len('THz.txt') # Total Number of Lines

    THz2eV = []
    with open('THz.txt', 'r') as f:
        lines = f.readlines()
        for line in lines:
            print(line.rstrip())
            words = line.split()
            THz2eV.append((int(words[0]), 0.001*float(words[-2])))  # Serial, convert to eV
    print(THz2eV)
    ####################################################

    #### Zero Point Energy (ZPE) ####
    #### += 0.5 * Summation of Frequency Energy [eV]  eV or meV
    
    #### Entropy Term (modified by JP)
    #                  x             
    # T*dS  =  kT [----------  -  ln(1 - exp(-x)) ]
    #              [exp(x)-1]
    # 
    #  x  = energy / kT = beta*hbar*omega
    #
    # <E(T)> = ZPE + kT ( x/[exp(x) - 1] )
    #
    # F_vib = <E(t)> -T<S(T)> = ZPE + kT *  ln(1 - exp(-x))
    #
    # K. Reuter et al. PRB (2001); S.-J. Woo et al. PRL (2013)
    ######################
    #for i in range(len(THz2eV)):
    ### Normal mode = 3 * Natom
    if natom:
        dof = natom * 3
    else:
        dof = len(THz2eV)

    kB = 0.0000861733576020577 # eV K-1
    T = 298.15 # Temeprature, K
    kT = kB * T

    Evib = 0.0      # Evib = E(T)
    TS = 0.0
    TS1 = 0.0
    TS2 = 0.0
    zpe = 0.0       # double of zpe
    #for i in range(len(THz2eV)):
    for i in range(dof):
        hw = THz2eV[i][1]       # eV
        zpe += hw               # sum of ZPE
        x = hw / kT
        den0 = math.exp(x) -1
        Evib += hw * (0.5 + 1./den0)

        ts1 = x / den0
        ts2 = -math.log(1. - math.exp(-x))
        ts = ts1 + ts2 
        TS += kT*ts
        TS1 += kT * ts1
        TS2 += kT * ts2
        #print(f"ts1 {ts1*kT:10.4f} ts2 {ts2*kT:10.4f}")
    zpe *= 0.5
    Fvib = Evib - TS    # == zpe - TS2
    #### Print Output ###
    print('----------------------------------------------')
    print(f'{"Zero Point Energy (ZPE)":>25} : {zpe:7.3f} eV (1)')
    print(f'{"Fvib 2nd logterm":>25} : {TS2:7.3f} eV (2)')
    print(f'{"E_vib (<E(T)>)":>25} : {Evib:7.3f} eV (3)')
    #print(f'{"Delta":>25} : {Evib-zpe:7.3f} eV')
    print(f'{"Entropy Energy (T<S>)":>25} : {TS:7.3f} eV (4)')
    #print(f'{"Entropy terms (TS1, TS2)":>25} : {TS1:7.3f} {TS2:7.3f} eV')
    print(f'{"Helmholtz free E (Fvib)":>25} : {Fvib:7.3f} eV', '(1)-(2) or (3)-(4)')
    print('----------------------------------------------')
    print(f"Don't forget to remember: Degree of Freedom == 3 *{natom} Natom ")
    with open("zpe.dat","w") as f:
        f.write('----------------------------------------------\n')
        f.write('Zero Point Energy (ZPE) : %15.9f  eV\n' % zpe)
        f.write('Entropy Energy (TS)     : %15.9f  eV\n' % TS)
        f.write('Helmholtz Energy (Fvib) : %15.9f  eV\n' % Fvib)
        f.write('----------------------------------------------\n')
    
    os.chdir(cwd)
    return 0

def main():

    parser = argparse.ArgumentParser(description='for ZPE & TS calculation to energy profile (HER, ORR, etc.)')
    parser.add_argument('d', nargs='?', default=os.getcwd(), help='input dir|OUTCAR which have OUTCAR w. ZPE calculation ')
    parser.add_argument('-na', '--natom', type=int, help='THz is redundant in OUTCAR so fix the number of tested atoms')
    args = parser.parse_args()

    cal_zpe_ts(args.d, args.natom)
    return 0

if __name__ == '__main__':
    main()



