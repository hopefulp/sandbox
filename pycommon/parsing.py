import re
import sys
from common import whereami
'''
def is_int_el
def is_int      
'''
### list

def find_letter(ch, lst):
    return any(ch in word for word in lst)

def find_letter2d(ch, lst2d):
    return any(ch in word for lst1d in lst2d for word in lst1d)

def convert_2lst_2Dlist(li, lshape):
    l2d = []
    for i in range(len(lshape)):
        lst = []
        while len(lst) < lshape[i]:
            print(f"{len(lst)} < {lshape[i]}")
            ### type of j: 1, 9, 11-15, 
            atom_0 = li.pop(0)
            if '-' in atom_0:
                ele = list(map(int, re.split('-', atom_0)))
                jlst = range(ele[0], ele[1]+1)
                lst.extend(jlst)
            else:
                lst.append(int(atom_0))
        l2d.append(lst)
    print(f"{whereami():>15}(): {l2d}")
    return l2d
    


def is_int_el(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def is_int(s):
    """ return number of integers if true otherwise False """

    _str = s.split()
    if _str:
        for ss in _str:
            if not is_int_el(ss):
                return False
        if len(_str) >= 1:
            return len(_str)
    else:
        print("null string")
        exit(1)

def is_there_char(st):
    """ return True if there is char """
    for char in st:
        if char.isalpha():
            return True
    return False
### string 
def str_decom(st, delimiter='.', index=0):
    lstr = st.split(delimiter)
    if len(lstr) < index+1:
        print("Error: %s has more than 1 dot in file name" % st)
        exit(1)
    else:
        return lstr[index]

def str_replace(st, old_str, new_str):
    newst = st
    return newst


def convert_s2l(st):
    list0=[]
    list0[:0] = st
    return list0

def get_atomlist4str(st):
    '''
    used in $SB/pyqchem/qcout_nbo.py
    '''
    list0=[]
    for ch in st:
        if ch.isupper():
            list0.append(ch)
        else:
            prech = list0.pop()
            newch = prech + ch
            list0.append(newch)
    return list0
    
