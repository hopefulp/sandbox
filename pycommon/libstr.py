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

  

   

