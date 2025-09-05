import sys, re, os
from parsing import anylower, f_parsing #, print_list2col
import numpy as np
from common import whereami

'''
pot_dict: list
    key word
    margin before potential
    Lx  False only extract pot (1 col) True x-pot (2 col)
'''
pot_dict = {'potcar': {'regi': 'local part', 'regf': 'lowercase', 'nskip': 1, 'x_coord': False }
            , 'siesta':  {'regi': 'Vna',     'regf': 'lowercase', 'nskip': 1, 'x_coord': True}
            }

def parse_potcar(pot, kw='local', opt=None):
    '''
    extract local part
    '''

    if pot_dict[kw]['x_coord']: Lx_coord = True     # extract 2 col       
    else:                       Lx_coord = False    # extract 1 col for only pot

    ### use f_parsing()
    lines, skip_lines = f_parsing (pot, pot_dict[kw]['regi'], pot_dict[kw]['regf'], nskip=pot_dict[kw]['nskip'])

    ### deal with pot_head
    if kw == 'local':
        head = float(skip_lines[0].strip())
    else:
        head = None

    ### check number of elements in a line
    eles = lines[0].strip().split()
    Lflatten = False
    if 2 < len(eles):           # if not x pot, 2 colums, all x or all pot -> flatten
        Lflatten = True

    ### parsing POTCAR
    pots = []

    for line in lines:
        ### 2D data converts into 1D
        if not Lx_coord and Lflatten:
            float_list = [float(item) for item in line.strip().split()]
            pots.extend(float_list)    # convert to 1D list
        else:
            potxy = line.strip().split()
            pots.append(potxy)
    #print(f"pots {pots}")
    #print(f"pots shape {np.array(pots).shape} in {whereami()}")
        
    return pots, head

        