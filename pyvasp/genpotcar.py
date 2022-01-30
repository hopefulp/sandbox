#!/home/joonho/anaconda3/bin/python

# Min Jong Noh
# Update : 2020 / 04 / 12

'''
This script for POTCAR generation.
POSCAR with element name at line 1 or 6 is required.
Potential directory should be located. (listed below)
'''

import os, sys, argparse

# Potential path
home = os.environ['HOME']
sb=home + "/sandbox/"
LDA_path  = sb + "pyvasp/vasppot/1.POTPAW.LDA.54.RECOMMEND"
#PBE_path  = "/TGM/Apps/VASP/POTCAR/2.POTPAW.PBE.54.RECOMMEND"
PBE_path  = sb + "pyvasp/vasppot/2.POTPAW.PBE.54.RECOMMEND"
PW91_path = sb + "pyvasp/vasppot/5.OUTDATED.POTPAW.GGA"

sv_list = ["Li", "K", "Rb", "Cs", "Sr", "Ba", "Sc", "Y", "Zr"]
pv_list = ["Na", "Ca", "Ti", "V", "Cr", "Mn", "Nb", "Mo", "Tc", "Hf", "Ta", "W", "Os"]
d_list  = ["Ga", "Ge", "In", "Sn", "Tl", "Pb", "Bi", "Po", "At"]

def fileread(fname):
    lineinfo = []
    wordinfo = []
    with open(fname) as f:
        for i, l in enumerate(f):
            line = l
            word = line.split()
            lineinfo.append(line)
            wordinfo.append(word)
    return lineinfo, wordinfo


def get_potcar(pot_type, hpp_list=None):

    error_check = 0
    if pot_type == 'lda':
        pot_path = LDA_path
    elif pot_type =='pbe':
        pot_path = PBE_path
    elif pot_type =='pw91':
        pot_path = PW91_path
    else:
        error_check = 1
        print("Warning, you should select among (lda, pbe, pw91)")

    # Read POSCAR
    pos_line, pos_word = fileread('POSCAR')

    # Listing up the element including VASP 4.X, 5.X ver.
    element = pos_word[5]
    element_check = element[0]

    try:
        int(element_check)
        element = pos_word[0]
    except ValueError:
        pass

    # Recommended VASP POTCAR

    element_refine = []

    for name in element:
        if name != 'H' or not hpp_list:
            if name in sv_list:
                element_refine.append(name+'_sv')
            elif name in pv_list:
                element_refine.append(name+'_pv')
            elif name in d_list:
                element_refine.append(name+'_d')
            else:
                element_refine.append(name)
        elif name == 'H' and hpp_list:
            hpp = hpp_list.pop(0)
            if hpp[0] == '1':
                hpp_dir = hpp[0] + '.' + hpp[1:]
            else:
                hpp_dir = '.' + hpp[:]
            element_refine.append(name+hpp_dir)
            print(f"find H directory {name+hpp_dir}")
                

    # Generation VASP POTCAR
    cmd = 'cat'
    if error_check == 1:
       print('Check what you are doing or to set a right POTCAR path')
    else:
       os.system('rm -rf POTCAR')
       for element in element_refine:
           addcmd = ' ' + '%s/%s/POTCAR' % (pot_path, element)
           cmd = cmd + addcmd 
       cmd = cmd + ' > POTCAR'
       os.system('%s' % cmd)
       print('%s POTCAR is sucessfully generated' % pot_type)

def main():
    parser = argparse.ArgumentParser(description='make POTCAR with pp flavor')
    parser.add_argument('-pp', '--pseudop', choices=['lda', 'pbe', 'pw91'], help='select pseudo potential')
    parser.add_argument('-hpp', '--pseudoH', nargs='*', help='pseudoH value')
    args = parser.parse_args()

    get_potcar(args.pseudop, args.pseudoH)

    return 0

if __name__ == '__main__':
    main()
