import argparse
import re
import os
import sys
import inspect
import numpy as np
import weakref
import glob
#from libstr import get_char_series_by_two
#from varname import nameof

"""
    version dependency included
    common functions gathering
    Class:: MyClass
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

class MyClass_obj(dict):
    instances = []                                              # class variable gathers object weak reference
    def __init__(self, name=None):
        self.__class__.instances.append(weakref.proxy(self))    # obj.__class__ : class MyClass, instance[] includes "self_object"
        self.name = name                                        # obj.name receives argument in a=MyClass(name)

class MyClass(MyClass_obj):
    pass

### is this working?
class MyClass_str(dict):
    instances = []      # gathers string
    def __init__(self, name=None):
        self.__class__.instances.append(name)                   # this is name, can't be used for obj.name 
        self.name = name

class CaseInsensitiveKey(object):
    def __init__(self, key):
        self.key = key
    def __hash__(self):
        return hash(self.key.lower())
    def __eq__(self, other):
        return self.key.lower() == other.key.lower()
    def __str__(self):
        return self.key

class CaseInsensitiveDict(dict):
    def __setitem__(self, key, value):
        key = CaseInsensitiveKey(key)
        super(CaseInsensitiveDict, self).__setitem__(key, value)
    def __getitem__(self, key):
        key = CaseInsensitiveKey(key)
        return super(CaseInsensitiveDict, self).__getitem__(key)


'''
    compare key 
    dic original dictionary
    key to find a key in dic
def compareKey(dic, key):
    keys = dic.keys()
    for k in keys:
        exp = re.compile(k, re.I)

'''

def get_digits_4str(word, string):
    ''' obtain digits after key-word '''
    if word in string:
        substr = re.split(word,string)[1]
        mat = re.match('\d+', substr)
        digits = mat.group()
    else:
        digits = None
    return digits


def dir_classify_n(lsorted, class_instance, class_dict,Lwrite=1):
    """
    classify files in a directory
    """
    #print(classobj_dict_key)
    luse=[]
    ukeys=[]        # used keys
    lsuff=['py','sh','csh','pl','']                 # '' include binary file without extension
    for f in class_dict.__dict__.keys():             # as for keys == py_fname(.py)
        Ltag = False

        for suf in lsuff:
            if suf:
                fsuf = f + '.' + suf
            ### include binary file without extenstion
            else:
                fsuf = f
            if fsuf in lsorted:                      # scan for keys
                luse.append(fsuf)
                lsorted.remove(fsuf)                 # remove key.py
                ukeys.append(f)                     # return key
                Ltag = True
                continue
        if Ltag ==  True:
            continue

    

        #print(f" in dir_classify():   {f}")      # to print all the not-selected files
    ### classify modules used
    if luse:
        CLASS_instance = class_instance.upper()
        print("  {:<10}::".format(CLASS_instance))
        if Lwrite:
            for f in luse:
                print(f"    {f}     ")
    return ukeys

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
def get_char_series_by_two(chars):
    '''
    chars should have two character
    '''
    ch1 = chars[0]
    ch2 = chars[-1]
    lchar = [chr(ch) for ch in range(ord(ch1), ord(ch2)+1)]
    return lchar
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

def whereami(rank=1):
    '''
    rank can be [0,1]
    '''
    return inspect.stack()[rank][3]

def yes_or_no(question):
    reply = str(input(question+' (y/n): ')).lower().strip()     # raw_input is renamed in v3.6
    print(reply)
    if re.match('y', reply):
        return True
    else:
        return False

def list2str(li, delimit=None):
    '''
    list to string
        delimit delimiter between list elements, default=''
    '''

    if delimit == None:
        st = "".join(str(x) for x in li)
    else:
        st = delimit.join(str(x) for x in li)
    return st

def list2dict(li):
    it = iter(li)
    dic = dict(zip(it,it))
    return dic

def print_list(li):
    st = list2str(li)
    print(st)
    return st

def get_answers(question):
    reply = str(input(question)).strip()
    #print reply
    return reply


def get_files_pattern(m_type, pattern,dirname, Ldir=False, Linverse=False, Lparents=None):
    #print(m_type)
    Lshow = False
    if m_type == 'p':
        f_list = get_files_prefix(pattern, dirname, Lshow, Ldir=Ldir, Lparents=Lparents)
    elif m_type == 's':
        f_list = get_files_suffix(pattern, dirname, Lshow, Ldir=Ldir, Linverse=Linverse, Lparents=Lparents)
    elif m_type == 'm':
        f_list = get_files_match(pattern, dirname, Lshow, Ldir=Ldir, Lparents=Lparents)
    else:
        print("No matching for file extraction")
    return f_list


def get_files_patterns(m_type, pattern, wdir, Ldir=False, Linverse=False, Lparents=None):
    """
    receive list of prefix & directory
    prefixes: list of prefix
    Lshow   to check list inside function
    """
    Lshow = False
    matched_files=[]
    ### Codes for prefix
    dir_files = os.listdir(wdir)
    i=0
    for fname in dir_files:
        for patt in pattern:
            if m_type == 'p':
                if re.match(patt, fname):
                    if Ldir or not os.path.isdir(fname):
                        matched_files.append(fname)
                    #print (patt, fname)
            ### for suffix
            elif m_type == 's':
                if fname.endswith(patt):
                    #if not Linverse:
                    if Ldir or not os.path.isdir(fname):
                        matched_files.append(fname)
                    ### included parents and directories
                    if Lparents:
                        #relative_files = get_relatives_suff(fname, dir_files)
                        fnlist = re.split('\.', fname)
                        if os.path.isdir(fnlist[0]):
                            rel_dir = fnlist[0]
                            print(f"{i:02d}: relative files {rel_dir}")
                            matched_files.append(rel_dir)
                            i += 1
            ### for search
            elif m_type == 'm':
                if re.search(patt, fname):
                ### if it is dir skip
                    if Ldir or not os.path.isdir(fname):
                        matched_files.append(fname)
            if Lshow:
                print(f"detect {fname}") # in {match} {matches}")
                
            #elif Linverse:
            #    if not os.path.isdir(fname):
            #        matched_files.append(fname)
    return matched_files
                
def get_files_type(ftype, dirname):
    matched_files=[]
    print(ftype)
    for fname in os.listdir(dirname):
        if fname.endswith(ftype):
            matched_files.append(fname)
    print(matched_files)
    return matched_files            

### deprecate get_dirs_prefix
def get_dirfiles(wdir, prefix=None, excludes=None, Lshow=True, Ldir=True):
    """
        receive list of prefix & directory
        prefixes: list of prefix
    """
    pwd = os.getcwd()
    matched_dirs=[]
    print(f'{__name__}', os.listdir(wdir))
    for fname in os.listdir(wdir):
        # re.match finds only prefix
        #print(f'1: {fname}')
        if os.path.isdir(f'{pwd}/{wdir}/{fname}'):      # need path for fname
            print(f'2: {fname}')
            if prefix:
                if not re.match(prefix, fname):
                    break       # to not list
                if excludes:
                    for ex in excludes:
                        if re.search(ex, fname):
                            break
            matched_dirs.append(fname)
            print (fname)

    return matched_dirs

def get_files_prefix(prefixes, dirname, Lshow=None, Ldir=None, sub_match=None):
    """
        receive list of prefix & directory
        prefixes: list of prefix
            sub_match   extract interval between a-g
        Ldir    T/F to include dir
    """
    matched_files=[]
    ### if sub_match, change prifixes list
    if sub_match:   # get two alphabet
        chars = get_char_series_by_two(sub_match)
        new_prefixes = []
        for pref in prefixes:
            for ch in chars:
                new_prefixes.append(pref + ch)
        prefixes = new_prefixes 

    for pref in prefixes:
        print(f"prefix: {pref} in {whereami()} of module {__name__}")
        for fname in os.listdir(dirname):
            # re.match finds only prefix
            if re.match(pref, fname):
                if not Ldir and os.path.isdir(fname):
                    continue
                matched_files.append(fname)
                if Lshow: print (pref, fname)
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
    return fname.split('/')[-1].split('.')[0]
fname_root = f_root
fname_ext  = f_ext

def fname_decom(fname):
    fname_parts = fname.split('.')
    if len(fname_parts) < 2:
        print("Error: %s has not dot in file name" % fname)
        exit(1)
    suff = fname_parts.pop()
    fn = '_'.join(fname_parts)
    return fn, suff        

def dirfname_decom(dirfname):
    dfname = dirfname.split('/')
    if len(dfname) != 2:
        print("Error: %s has not dot in file name" % fname)
        exit(1)
    dirname = dfname[0]
    fname_parts = fname_decom(dfname[1])
    suff = fname_parts.pop()
    fn = '_'.join(fname_parts)
    return dirname, fname_parts[0], suff, fn       

fname_parsing = fname_decom

def find_file(dname, f_root):
    files = os.listdir(dname)
    for f in files:
        if f.endswith(f_root):
            return f
        elif f.startswith(f_root):
            return f
    return None

def get_files_suffix_list(suffixes, flist, Lshow=False, Ldir=False):
    """
        receive list of suffixes & files 
        suffixes : list of suffix
    """
    matched_files=[]
    dirs=[]
    files=[]
    for fname in flist:
        if os.path.isdir(fname):
            dirs.append(fname)
        else:
            files.append(fname)
    for suff in suffixes:
        for fname in files:
            #print(f" {suff} in {fname} ?")
            if fname.endswith(suff):
                matched_files.append(fname)
    matched_files.extend(dirs)                
    return matched_files    

def get_files_suffix(suffixes, dname, Lshow=False, Ldir=False, Linverse=False, Lparents=False):
    """
        receive list of suffixes & directory
        suffixes : list of suffix
    """
    matched_files=[]
    for suff in suffixes:
        for fname in os.listdir(dname):
            ### if matched
            if fname.endswith(suff):
                ### if Ldir, include dir || if not dir, include
                if not Linverse:
                    if Ldir or not os.path.isdir(fname):
                        matched_files.append(fname)
            elif Linverse:
                if not os.path.isdir(fname):
                    matched_files.append(fname)
    return matched_files                

def get_files_match(matches, dname, Lshow, Ldir=False):
    """
        receive list of patterns for match
        matches: list of match
    """
    matched_files=[]
    ### two for-loop is converted
    for fname in os.listdir(dname):
        #print(f"{fname}")
        for match in matches:
            if re.search(match, fname):
                ### if it is dir skip
                if os.path.isdir(fname) and not Ldir:
                        continue
                if Lshow:
                    print(f"detect {fname}") # in {match} {matches}")
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

def search_dirs(dir_pre, fname):
    dirs =  glob.glob(f"{dir_pre}*")
    list_dir=[]
    #print(fname, dirs)
    for d in dirs:
        if os.path.isfile(f"{d}/{fname}"):
            list_dir.append(d)
    return list_dir            
        

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

def fname_index(pre, suff, cwd):
    files = os.listdir(cwd)
    for f in files:
        if re.match(pre, f) and re.search(suff, f):
            pass
    return 0

def main():
    parser = argparse.ArgumentParser(description='test functions')
    parser.add_argument('-d', '--dirs', help='dirname')
    parser.add_argument('-f', '--funcname', help='function name')
    args = parser.parse_args()

    get_dirs_prefix(args.dirs)

if __name__ == '__main__':
    main()
