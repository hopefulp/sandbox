"""
    common functions gathering
    yes_or_no(question) : to obtain y/n for next command
    get_files_prefix(list_prefix, directory) : returns list
    get_files_suffix(list_suffix, directory) : returns list
    fname_decom(fname): returns prefix, extension
"""

import re
import os
import inspect

def whereami():
    return inspect.stack()[1][3]

def yes_or_no(question):
    reply = str(raw_input(question+' (y/n): ')).lower().strip()
    if re.match('y', reply):
        return True
    else:
        return False

def get_answers(question):
    reply = str(raw_input(question)).strip()
    #print reply
    return reply


def get_files_pattern(m_type, pattern, dir):
    if m_type == 'p':
        f_list = get_files_prefix(pattern, dir)
    elif m_type == 's':
        f_list = get_files_suffix(pattern, dir)
    elif m_type == 'm':
        f_list = get_files_match(pattern, dir)
    return f_list

def get_files_prefix(prefixes, dir):
    """
        receive list of prefix & directory
        prefixes: list of prefix
    """
    matched_files=[]
    for pref in prefixes:
        for file in os.listdir(dir):
            # re.match finds only prefix
            if re.match(pref, file) and not os.path.isdir(file):
                matched_files.append(file)
                #print pref, file
    return matched_files

def get_files_suffix(suffixes, dir):
    """
        receive list of suffixes & directory
        suffixes : list of suffix
    """
    matched_files=[]
    for suff in suffixes:
        for file in os.listdir(dir):
            if file.endswith(suff) and not os.path.isdir(file):
                matched_files.append(file)
    return matched_files                

def get_files_match(matches, dir):
    """
        (not checked) receive list of patterns for match
        matches: list of match
    """
    matched_files=[]
    for match in matches:
        for file in os.listdir(dir):
            if re.search(match, file) and not os.path.isdir(file):
                matched_files.append(file)
    return matched_files

def get_files_exclude(matches, dir):
    """
        if not matched,
        return files
    """
    all_files=os.listdir(dir)
    save_list=[]
    imatch=0
    for match in matches:
        for file in all_files:
            ### exclude dir once
            if imatch==0:
                if os.path.isdir(file):
                    save_list.append(file)
                    #print file
                    continue
            if re.search(match, file):
                save_list.append(file)
                #print file
        imatch+=1                
    for file in save_list:
        #print file
        all_files.remove(file)
    return all_files

def fname_decom(fname):
    lname = fname.split('.')
    if not len(lname) == 2:
        print("Error: %s has more than 1 dot in file name" % fname)
        exit(1)
    return lname[0], lname[1]        
                
