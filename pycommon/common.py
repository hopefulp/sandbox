import re
import os
import sys
import inspect
import numpy as np

"""
    version dependency included
    common functions gathering
    yes_or_no(question) : to obtain y/n for next command
    FILES
        find_file(dir,kw): search dir & find fname with kw in front or end
        dir_files(): returns files with execution|module
        get_files_prefix(list_prefix, directory) : returns list
        get_files_suffix(list_suffix, directory) : returns list
        fname_decom(fname): returns prefix, extension
        f_ext(fname): get extension
        f_root(fname): get root
    whereami(): returns function name
"""

class MyClass(dict):
    pass


def dir_classify(lsorted, classobj_dict_key,classobj_dict,Lwrite=1):
    """
    classify files in a directory
    """
    #print(classobj_dict_key)
    c_obj = classobj_dict[classobj_dict_key]
    luse=[]
    ukeys=[]        # used keys
    for f in c_obj.__dict__.keys():             # as for keys == py_fname(.py)
        fpy=f+".py"
        if fpy in lsorted:                      # scan for keys
            luse.append(fpy)
            lsorted.remove(fpy)                 # remove key.py
            ukeys.append(f)                     # return key
            continue
        #print(f" in dir_classify():   {f}")      # to print all the not-selected files
    ### classify modules used
    if luse:
        print("  {:<10}::".format(classobj_dict_key+" used"))
        if Lwrite:
            for f in luse:
                print(f"    {f}     ")
    return ukeys

def classify_dirs(lsorted, classobj_dict_key,classobj_dict):
    """ classify directories """
    c_obj = classobj_dict[classobj_dict_key]
    luse=[]
    ukeys=[]

    #print(c_obj.__dict__.keys())
    #print(lsorted)
    for d in c_obj.__dict__.keys():
        if d in lsorted:
            luse.append(d)
            lsorted.remove(d)       # lsorted affect the calling module
            #ukeys.append(d)
    ### classify modules used
    if luse:
        print(f"  {classobj_dict_key} -- Dirs ::")
        for d in luse:
            print(f"    {d}")
    return luse

def whereami():
    return inspect.stack()[1][3]

def yes_or_no(question):
    reply = str(input(question+' (y/n): ')).lower().strip()     # raw_input is renamed in v3.6
    print(reply)
    if re.match('y', reply):
        return True
    else:
        return False

def list2str(li):
    st = "".join(str(x) for x in li)
    return st

def print_list(li):
    st = list2str(li)
    print(st)
    return st

def get_answers(question):
    reply = str(input(question)).strip()
    #print reply
    return reply


def get_files_pattern(m_type, pattern,dirname):
    #print(m_type)
    if m_type == 'p':
        f_list = get_files_prefix(pattern, dirname)
    elif m_type == 's':
        f_list = get_files_suffix(pattern, dirname)
    elif m_type == 'm':
        f_list = get_files_match(pattern, dirname)
    else:
        print("No matching for file extraction")
    return f_list

def get_files_type(ftype, dirname):
    matched_files=[]
    print(ftype)
    for fname in os.listdir(dirname):
        if fname.endswith(ftype):
            matched_files.append(fname)
    print(matched_files)
    return matched_files            

def get_files_prefix(prefixes, dirname):
    """
        receive list of prefix & directory
        prefixes: list of prefix
    """
    matched_files=[]
    for pref in prefixes:
        for fname in os.listdir(dirma,e):
            # re.match finds only prefix
            if re.match(pref, fname) and not os.path.isdir(fname):
                matched_files.append(fname)
                #print pref, file
    return matched_files

def dir_files(dir_):
    '''
    separate files into executable and module in that directory
    '''
    executable=[]
    mod=[]
    for f in os.listdir(dir_):
        if re.match('\.', f) or os.path.isdir(dir_+'/'+f) or f_ext(f) == 'pyc' or f_ext(f) == 'bak' :
            pass
        elif os.access(dir_+'/'+f, os.X_OK):
            executable.append(f)
        else:
            mod.append(f)
    return executable, mod            

def dir_all(dir_):
    '''
    all the directories and files
    '''
    dirs=[]
    executable=[]
    mod=[]
    d_link_indices=[]
    for f in os.listdir(dir_):
        if re.match('\.', f) or f_ext(f) == 'pyc' or f_ext(f) == 'bak' :
            pass
        elif os.path.isdir(os.path.join(dir_, f)):
            dirs.append(f)
            if os.path.islink(os.path.join(dir_, f)):
                d_link_indices.append(dirs.index(f))
        elif os.access(os.path.join(dir_, f), os.X_OK):
            executable.append(f)
        else:
            mod.append(f)
    return executable, mod, dirs, d_link_indices

def f_ext(fname):
    return fname.split('.')[-1]
def f_root(fname):
    return fname.split('.')[0]
fname_root = f_root
fname_ext  = f_ext

def fname_decom(fname):
    lname = fname.split('.')
    if not len(lname) == 2:
        print("Error: %s has more than 1 dot in file name" % fname)
        exit(1)
    return lname[0], lname[1]        

fname_parsing = fname_decom

def find_file(dname, f_root):
    files = os.listdir(dname)
    for f in files:
        if f.endswith(f_root):
            return f
        elif f.startswith(f_root):
            return f
    return None

def get_files_suffix(suffixes, dname):
    """
        receive list of suffixes & directory
        suffixes : list of suffix
    """
    matched_files=[]
    for suff in suffixes:
        for fname in os.listdir(dname):
            if fname.endswith(suff) and not os.path.isdir(fname):
                matched_files.append(fname)
    return matched_files                

def get_files_match(matches, dname, Lshow):
    """
        receive list of patterns for match
        matches: list of match
    """
    matched_files=[]
    for match in matches:
        for fname in os.listdir(dname):
            if re.search(match, fname) and not os.path.isdir(fname):
                if Lshow:
                    print("detect {fname}")
                matched_files.append(fname)
    return matched_files

def get_files_exclude(matches, dname):
    """
        if not matched,
        return files
    """
    all_files=os.listdir(dname)
    save_list=[]
    imatch=0
    for match in matches:
        for fname in all_files:
            ### exclude dir once
            if imatch==0:
                if os.path.isdir(fname):
                    save_list.append(fname)
                    #print file
                    continue
            if re.search(match, fname):
                save_list.append(fname)
                #print file
        imatch+=1                
    for fname in save_list:
        #print file
        all_files.remove(fname)
    return all_files



def expand_dim_str(lstring):
    """ expand 1D string to 2D string """
    a = np.array(lstring)
    dim = len(a.shape)
    new_2d=[]
    if dim == 1:
        ### ele is string
        for ele in a:
            new_list=ele.split()
            new_2d.append(new_list)
    return new_2d

def f_parsing(fname, sep='num'):
    with open(fname, 'r') as f:
        lines = f.readlines()
        i=0
        y=[]
        for line in lines:
            if re.search('[a-zA-Z]', line):
                pass
            else:
                i += 1
                line = line.strip()
                ele = line.split()
                if i == 1:
                    n = len(ele)
                    for j in range(n):
                        y.append([])         # declare list as many as the number of column
                for j in range(n):
                    y[j].append(float(ele[j]))
    for j in range(n):
        if j > 0:
            if len(y[j-1]) != len(y[j]):
                print("in zip file to column, Error: length of colums are different")
                sys.exit(1)
    return y
