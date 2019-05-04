import numpy as np
import re
import os
from common import *
from qcout_ini import *
from my_print import *
import heapq
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
            print ("Error:: No atom filtering type in", whereami() )
            exit(55)
 
def f_imo_basis(imoc_list, l_atoms, filter_tag, mo_id):
    """
    decide whether 1 level of dump block of 'imoc_list' as for the given mo_id is saved for draw
    imoc_list is saved for the atoms listed in atom_list
    """
    
    #### change element which corresponds to l_atoms ordering, to check its existence with high MO coefficients
    n_latoms = len(l_atoms)
    atom_tag=np.zeros(n_latoms)
    base_coeff = {}
    mcoeff = 0.0
    #print l_atoms, "in", whereami()
    if not imoc_list:
        return 0, base_coeff
    else:
        k=0
        ### Obtain atom_tag: which atom is included in imo one block
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
        
        ### atom existence was checked in atom_tag
        #   for SUB, atom_tag == [1,0,0,1] etc: 1st atom should exist
        #print atom_tag, "in", mo_id, "with", filter_tag, "in", whereami()
        if filter_tag=="ONE" :
            if atom_tag[0]==1:
                base_coeff = f_imo_basis_dic(imoc_list)
                selection= 1
            else: 
                selection= 0
        elif re.search("SUB", filter_tag):
            if atom_tag[0]==1: # Ni is neccessary 
                base_coeff = f_imo_basis_dic(imoc_list)
                selection= 1
            else:
                selection= 0
        elif re.search("SEL", filter_tag):
            if atom_tag[0]==1:
                base_coeff = f_imo_basis_dic(imoc_list)
                selection= 1
            else:
                selection= 0
        elif filter_tag=="ALL":             
                base_coeff = f_imo_basis_dic(imoc_list)
                selection= 1
        else:
            print ("Error:: No atom filtering type in", whereami() )
            exit(55)
    #print base_coeff
    new_dic = trim_coeff(base_coeff, n_latoms*Nbasis_show)
    if V_print_moc >=1: 
        print ("%5d in %s" % (mo_id, whereami() ))
        lprint_sorted(new_dic)
    return selection, new_dic

def trim_coeff(dic, n):
    """ trim dict to leave some maximum values """
    if V_print_moc >= 2:
        print (whereami())
        lprint_sorted(dic)
    if len(dic) > n:
        a = heapq.nlargest(n, dic.items(), lambda i: i[1])
        if V_print_moc >= 2: lprint(a)    
        new_dic = dict(a)
        return new_dic
    else:
        return dic

def f_imo_basis_dic(imoc_list):
    i=0
    bcoeff={}
    for line in imoc_list:
        if i == 0:
            if V_print_moc >= 2: print (line)
        else:
            if V_print_moc >= 2: print (line)
            b_coeff = line.split()
            coeff = abs(float(b_coeff[1]))
            if not b_coeff[0] in bcoeff.keys():
                bcoeff[b_coeff[0]] = coeff
            else:
                bcoeff[b_coeff[0]] += coeff
        i+=1
    if i==1:
        pass

    if V_print_moc >= 2: print ("bases_coeff dictionary:", bcoeff  , whereami())
    return bcoeff

def f_imo_dump1(moc_list, l_atoms, filter_tag, mo_id):
    """
    decide whether 1 level of dump block of 'moc_list' as for the given mo_id is saved for draw
    moc_list is saved for the atoms listed in atom_list
    """

    #### change element which corresponds to l_atoms ordering, to check its existence with high MO coefficients
    atom_tag=np.zeros(len(l_atoms))
    if not moc_list:
        return 0
    else:
        f_imo_print(moc_list)
        return 1

 
def f_imo_maxcoeff(moc_list):
    i=0
    max_coeff = 0
    for line in moc_list:
        if i == 0:
            print (line)
        else:
            basis_coeff = line.split()
            coeff = abs(float(basis_coeff[1]))
            if max_coeff < coeff:
                max_coeff = coeff
                max_basis = basis_coeff[0]
        if V_print >= 1: print  (line)
        i+=1
    if i==1:
        max_basis = 'NONE'
        max_coeff = 0
    print ("     max_basis: ",  max_basis, max_coeff, whereami())
    return max_basis, max_coeff


def f_imo_print(moc_list):
    for line in moc_list:
        if V_print_moc >= 1: print  (line)
    return 0


