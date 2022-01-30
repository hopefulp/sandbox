from __future__ import with_statement
import re
import os

"""
py module to extract lines from a file
"""

def extract_one_line(fname, kw, nline=1, opt=1):
    """
    if there is keyword, get lines
    call f(fname, kw): returns one line [nline=1], from next line [opt=1]
    """
    tag = "OFF"
    lines=[]        # extend to multiple lines
    try:
        with open(fname, 'r') as f:
            for line in f:
                if tag == "OFF":
                    if re.search(kw, line):
                        tag = "ON"
                # print next line
                else:
                    return line
    except IOError:
        print('error in %s' % __file__)
    return

def extract_values(fname, key):
    
    values=[]
    with open(fname, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if key in line:
                ele = line.strip().split()
                ind = ele.index(key)
                values.append(ele[ind+2])

    return values


def extract_col_qchem(fname, kw):
    """
    extract column list for atom name
    """
    tag = "OFF"
    aname = []
    with open(fname, 'r') as f:
        for line in f:
            if tag == "OFF":
                if re.search(kw, line):
                    i = 0
                    tag = "ON"
            elif tag == "ON":
                i += 1
                if i == 1: continue
                if re.match("\w", line):
                    line_ = line.split()
                    aname.append(line_[0])
                else:
                    return aname
    return                    
