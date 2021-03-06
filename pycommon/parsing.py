import re

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
