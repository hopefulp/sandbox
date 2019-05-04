import numpy as np
import re
import os
from common import *
from mplt_mo_ini import *
check_id=130
def f_imo_dump(moc_list, l_atoms, filter_tag, mo_id):
    """
    decide whether 1 level of dump block of 'moc_list' as for the given mo_id is saved for draw
    moc_list is saved for the atoms listed in atom_list
    """
    
    #### change element which corresponds to l_atoms ordering, to check its existence with high MO coefficients
    #print l_atoms, "in", whereami()
    atom_tag=np.zeros(len(l_atoms))
    if not moc_list:
        return 0
    else:
        k=0
        for atom in l_atoms:
            flag="OFF"
            j=0
            for moc_line in moc_list:
                #### skip first line of "MO level  mo-id  mo-ene"
                #print moc_line, "in", whereami()
                if j==0:
                    j+=1
                    continue
                #### in C-1-py 0.27933, we can find atom species in capital C,
                #### py, s, dx, f etc is lowercase
                #### do not use search for N as for "N Ni"
                fields=re.split("\s", moc_line)
                fields=[x for x in fields if x]
                t_atom_basis=re.split("-", fields[0])
                atom_id=t_atom_basis[0]
                #if re.search(atom, moc_line):
                #print atom, atom_name
                if atom == atom_id:
                    flag="ON"
                    #print atom, " in ", moc_line
                    break
            if flag=="ON":
                atom_tag[l_atoms.index(atom)]=1

        #### atom existence was checked in atom_tag
        #print atom_tag, "in", mo_id, "with", filter_tag, "in", whereami()
        if filter_tag=="ONE" :
            if atom_tag[0]==1:
                f_imo_print(moc_list)
                return 1
            else: 
                return 0
        elif re.search("SUB", filter_tag):
            if atom_tag[0]==1:
                f_imo_print(moc_list)
                return 1
            else:
                return 0
        elif re.search("SEL", filter_tag):
            if atom_tag[0]==1:
                f_imo_print(moc_list)
                return 1
            else:
                return 0
        elif filter_tag=="ALL":             
                f_imo_print(moc_list)
                return 1
        else:
            print "Error:: No atom filtering type in", whereami() 
            exit(55)
 

def f_imo_print(moc_list):
    for line in moc_list:
        if V_print_moc >= 1: print  line
    return 0


