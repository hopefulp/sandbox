#!/home/joonho/anaconda3/bin/python

### Update : 2019 / 03 / 08
# Min Jong Noh #
# starnmj@kaist.ac.kr

##################################################################################
# This is a script   #
##################################################################################

###
import os, sys, math
import argparse
###


###################################

#### Get a total number of lines in read-file ####

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i+1

def cal_zpe_ts(dir1, natom, infile):
    cwd = os.getcwd()
    fname = dir1 + '/' + infile
    if not os.path.isfile(fname):
        print(f"There is no {infile} in {dir1}")
        sys.exit(1)

    os.chdir(dir1)
    #### THz Data mining in OUTCAR ####
    os.system('grep THz OUTCAR > THz.txt')

    TN_lines = file_len('THz.txt') # Total Number of Lines

    THz2meV = []
    with open('THz.txt', 'r') as f:
        lines = f.readlines()
        for line in lines:
            print(line.rstrip())
            words = line.split()
            THz2meV.append((int(words[0]), 0.001*float(words[-2])))  # Serial, eV

    ####################################################


    #### Zero Point Energy (ZPE) ####
    #### += 0.5 * Summation of Frequency Energy [eV] ####

    ZPE = 0.0
    #for i in range(len(THz2meV)):
    ### Normal mode = 3 * Natom
    if natom:
        dof = natom * 3
    else:
        dof = len(THz2meV)
    for i in range(dof):
        E = 0.5 * THz2meV[i][1]
        ZPE += E
    ##################################

    #### Entropy Term ####
    #              x             
    # T*dS  =  ----------  -  ln(1 - exp(-x))
    #          [exp(x)-1]
    # 
    #  x  = energy / kT
    #
    ######################
    kB = 0.0000861733576020577 # eV K-1
    T = 298.15 # Temeprature, K
    kT = kB * T

    TS = 0.0
    #for i in range(len(THz2meV)):
    for i in range(dof):
        x = THz2meV[i][1] / kT
        v1 = x / (math.exp(x) - 1)
        v2 = 1 - math.exp(-x)
        E = v1 - math.log(v2)
        TS += kT*E

    #### Print Output ###
    print('----------------------------------------------')
    print('Zero Point Energy (ZPE) : %15.9f  eV' % ZPE)
    print('Entropy Energy (TS)     : %15.9f  eV' % TS)
    print('----------------------------------------------')
    print("Don't forget to remember: Degree of Freedom == 3 * Natom")
    with open("zpe.dat","w") as f:
        f.write('----------------------------------------------\n')
        f.write('Zero Point Energy (ZPE) : %15.9f  eV\n' % ZPE)
        f.write('Entropy Energy (TS)     : %15.9f  eV\n' % TS)
        f.write('----------------------------------------------\n')
    
    os.chdir(cwd)
    return 0

def main():

    parser = argparse.ArgumentParser(description='for ZPE & TS calculation to energy profile (HER, ORR, etc.)')
    parser.add_argument('d', nargs='?', default=os.getcwd(), help='input dir which have OUTCAR w. ZPE calculation')
    parser.add_argument('-na', '--natom', type=int, help='THz is redundant in OUTCAR so fix the number of tested atoms')
    parser.add_argument('-inf', '--infile', default='OUTCAR', help='OUTCAR')
    args = parser.parse_args()

    cal_zpe_ts(args.d, args.natom, args.infile)
    return 0

if __name__ == '__main__':
    main()



