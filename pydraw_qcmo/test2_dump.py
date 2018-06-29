import numpy as np
import re
import os
from common import *
from mplt_mo_ini import *
check_id=130
 
def f_imo_basis(imoc_list, l_atoms, filter_tag, mo_id):
    """
    decide whether 1 level of dump block of 'imoc_list' as for the given mo_id is saved for draw
    imoc_list is saved for the atoms listed in atom_list
    """
    
    #### change element which corresponds to l_atoms ordering, to check its existence with high MO coefficients
    atom_tag=np.zeros(len(l_atoms))
    mbasis = {}
    mcoeff = 0.0
    print l_atoms, "in", whereami()
    if not imoc_list:
        return 0, mbasis
    else:
        k=0
        for atom in l_atoms:
            flag="OFF"
            j=0
            for moc_line in imoc_list:
                #### skip first line of "MO level  mo-id  mo-ene"
                #print moc_line, "in", whereami()
                if j==0:
                    j+=1
                    continue
                #### in C-1-py 0.27933, we can find atom species in capital C,
                #### py, s, dx, f etc is lowercase
                #### do not use search for N as for "N Ni"
                fields=moc_line.split()
                fields=[x for x in fields if x]
                atom_type=re.split("-", fields[0])
                atom_name=atom_type[0]
                #if re.search(atom, moc_line):
                #print atom, atom_name
                if atom == atom_name:
                    flag="ON"
                    #print atom, " in ", moc_line
                    break
            if flag=="ON":
                atom_tag[l_atoms.index(atom)]=1

        #### atom existence was checked in atom_tag
        #print atom_tag, "in", mo_id, "with", filter_tag, "in", os.path.basename(__file__)
        if filter_tag=="ONE" :
            if atom_tag[0]==1:
                mbasis, mcoeff = f_imo_print(imoc_list)
                selection= 1
            else: 
                selection= 0
        elif re.search("SUB", filter_tag):
            if atom_tag[0]==1:
                mbasis, mcoeff = f_imo_print(imoc_list)
                selection= 1
            else:
                selection= 0
        elif re.search("SEL", filter_tag):
            if atom_tag[0]==1:
                mbasis, mcoeff = f_imo_print(imoc_list)
                selection= 1
            else:
                selection= 0
        elif filter_tag=="ALL":             
                mbasis, mcoeff = f_imo_print(imoc_list)
                selection= 1
        else:
            print "Error:: No atom filtering type in", whereami() 
            exit(55)
    #print selection, mbasis, mcoeff
    return selection, mbasis                


def f_imo_print(moc_list):
    for line in moc_list:
        if V_print_moc >= 1: print  line
    return 0


