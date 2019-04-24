import re
import os
import inspect

"""
    common functions gathering
    yes_or_no(question) : to obtain y/n for next command
    dir_files(): returns files with execution|module
    get_files_prefix(list_prefix, directory) : returns list
    get_files_suffix(list_suffix, directory) : returns list
    fname_decom(fname): returns prefix, extension
    whereami(): returns function name
"""

def whereami():
    return inspect.stack()[1][3]

def yes_or_no(question):
    #reply = str(raw_input(question+' (y/n): ')).lower().strip()
    reply = str(input(question+' (y/n): ')).lower().strip()     # raw_input is renamed in v3.6
    print(reply)
    if re.match('y', reply):
        return True
    else:
        return False

def get_answers(question):
    reply = str(raw_input(question)).strip()
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

def get_files_match(matches, dname):
    """
        (not checked) receive list of patterns for match
        matches: list of match
    """
    matched_files=[]
    for match in matches:
        for fname in os.listdir(dname):
            if re.search(match, fname) and not os.path.isdir(fname):
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

def fname_decom(fname):
    lname = fname.split('.')
    if not len(lname) == 2:
        print("Error: %s has more than 1 dot in file name" % fname)
        exit(1)
    return lname[0], lname[1]        
                
def fname_parsing(fname):
    lname = fname.split('.')
    if not len(lname) == 2:
        print("Error: %s has more than 1 dot in file name" % fname)
        exit(1)
    return lname[0], lname[1]        
