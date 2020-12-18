### string-related function

def convert_s2l(st):
    list0=[]
    list0[:0] = st
    return list0

def get_atomlist4str(st):
    list0=[]
    for ch in st:
        if ch.isupper():
            list0.append(ch)
        else:
            prech = list0.pop()
            newch = prech + ch
            list0.append(newch)
    return list0            
