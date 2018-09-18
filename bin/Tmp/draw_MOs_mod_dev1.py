import numpy as np
import re

def f_imo_dump(moc_list, l_atoms, filter_tag):
    if not moc_list:
        return 0
    else:
        i=0
        for atom in l_atoms:
            flag="OFF"
            j=0
            for moc_line in moc_list:
                #print "in f_imo_dump:: ", moc_line
                #### skip first line of "MO level  mo-id  mo-ene"
                if j==0:
                    j+=1
                    continue
                if re.search(atom, moc_line):
                    '''
                    if atom=="Ni":
                        if re.search("d", moc_line):
                            flag="ON"
                            break
                    else:
                    '''
                    flag="ON"
                    #print atom, " in ", moc_line
                    break
            if flag=="ON":
                i+=1
        if filter_tag=="ONE":
            if i>=1:
                f_imo_print(moc_list)
                return 1
            else: return 0
        else: 
            if i>=2:
                f_imo_print(moc_list)
                return 1
            else: return 0
   
def f_imo_print(moc_list):
    for line in moc_list:
        print line
    return 0

