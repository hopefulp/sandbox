#!/home/joonho/anaconda3/envs/py310/bin/python

import argparse
import os
import re
import sys

from libanalysis import defect_band_electron_count

def obtain_vbm(input):
    from pymatgen.io.vasp import Vasprun

    vr = Vasprun(input, parse_potcar_file=False)
    bs = vr.get_band_structure()
    print("VBM:", bs.get_vbm())
    print("CBM:", bs.get_cbm())
    print("Band gap:", bs.get_band_gap())



def main():
    '''
    job     vbm         pymatgen runs above python 3.10 and pymatgen.__version__
            band charge 
    '''
    parser = argparse.ArgumentParser(description='How to make a continuous job dir')
    ### job
    parser.add_argument('-j', '--job', default='band_charge', choices=['band_charge','vbm'], help='VASP extra jobs')
    parser.add_argument('-i', '--index', type=int, help='get integer')
    
    args = parser.parse_args()

    if args.job == 'vbm':
        obtain_vbm('vasprun.xml')
    elif args.job == 'band_charge':
        n = defect_band_electron_count("EIGENVAL", target_band=args.index) #
        print(f"Electrons in band {args.index} =", n)

    
if __name__ == '__main__':
    main()
