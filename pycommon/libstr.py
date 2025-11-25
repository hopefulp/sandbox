'''
string modification
imported from libfile, 

def isnumber(str) 
            returns True|False whether string is a number (float)
    obtain_sep(str)
            in  str is line
            return sep(arator) for delimiter

'''
from common import whereami

Lprint = 0
def isnumber(s):
    '''
    returns [True|False]
    '''
    try:
        float(s)
        return True
    except ValueError:
        return False  
    
def find_sep(s):
    '''
    search & return ',', '\t' in order, if not, returns None
    '''
    if "," in s:
        sep = ","
        if Lprint: print(f"delimiter is {sep}")
    elif "\t" in s:
        sep = "\t"
        if Lprint: print(f"delimiter is tab")
    else:
        sep = None   # let pandas auto-detect (whitespace)
        if Lprint: print(f"delimiter is None in {whereami()}")
    return sep

def li2str(li, delimit=None, Lline=False):
    '''
    list to string
        delimit delimiter between list elements, default=''
    '''

    if Lline:
        endline='\n'
    else:
        endline=""

    if delimit == None:
        st = "".join(str(x) for x in li) + endline
    else:
        st = delimit + delimit.join(str(x) for x in li) + endline
    return st
'''
def li2line(li, delimit=None):
    
    #new_line = " ".join(latoms) + "\n"
    new_line = " ".join(str(x) for x in li) + "\n"
    return new_line
'''
def list2dict(li):
    it = iter(li)
    dic = dict(zip(it,it))
    return dic

def print_list(li):
    st = list2str(li)
    print(st)
    return st



   

