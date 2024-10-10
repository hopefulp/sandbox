import re
#import sys
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

def convert_2lst2D(li, lshape):
    #print(f"{li} {lshape}")
    l2d = []
    tmplst = []
    for shap in lshape:
        lst = []
        while len(lst) < shap:
            #print(f"{len(lst)} < {lshape[i]} in module {__name__}")
            ### type of j: 1, 9, 11-15,  sometimes 0-100
            ### if not tmplost get from li
            if not tmplst and li:
                list_ele = li.pop(0)   # --> error when lack of li elements
                ### if -1, save -1 for Tdos: modified for '-1' of Tdos
                if list_ele == '-1':
                    tmplst=[-1]
                elif '-' in list_ele:
                    ele = list(map(int, re.split('-', list_ele)))
                    tmplst = list(range(ele[0], ele[1]+1))
                elif list_ele.isdigit():
                    tmplst.append(int(list_ele))
                else:
                    print("Error:: input type of atom list")
            #print(f"a: tmplist {tmplst}")
            ### tmplst exists
            nreq = shap - len(lst)
            if len(tmplst) <= nreq:
                lst.extend(tmplst)
                tmplst=[]
            else:
                lst.extend(tmplst[0:nreq])    # use slicing and del
                del tmplst[:nreq]
            #print(f"b: tmplist {tmplst}")
        l2d.append(lst)
    print(f"{whereami():>15}(): {l2d}")
    return l2d
    
### word

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

def startnum(word):
    '''return the index of the first number in word'''
    m = re.search("\d", word)
    #print(f'm.start() returns {m.start()}')
    if m:
        return m.start()
    else:
        return -1


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
    
