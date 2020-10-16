#!/home/joonho/anaconda3/bin/python

import argparse
import re
import importlib

def jobs(mod_comm, att, subkey):
    #subkeys = mod_comm.__dict__[att].__dict__.keys()
    print(att)
    subkeys = list(mod_comm.__dict__[att].__dict__) # list of dict returns keys
    print('---------------------------------------------')
    print(f"key:: {att}, sub_keys:: {subkeys}")
    sel_subkeys=[]
    #if isinstance(mod_comm.__dict__[att].__dict__[keyss[0]],dict):
    ### use select next key here
    if subkey:
        #subkey = subkey.upper()
        sel_subkeys.append(subkey)
    else:
        sel_subkeys.extend(subkeys[:])
    for key in sel_subkeys:
        #print(key)          #1st att(server), 2nd sge 
        if isinstance(mod_comm.__dict__[att].__dict__[key], dict):   # value of 'sge' is dict
        #print("here is True")
            for k in mod_comm.__dict__[att].__dict__[key].__dict__.keys():
                #print(f" k is {k}")
                print(mod_comm.__dict__[att].__dict__[key].__dict__[k])
            #sub_keys.append(key)
        else:
            print(mod_comm.__dict__[att].__dict__[key])
    if not subkey:
        print(f"Use -k for subkey in {subkeys}")

    return 0

def main():
    parser = argparse.ArgumentParser(description="shows dictionary for all: work, system, package  ")
    #parser.add_argument('-m', '--mod', default='sys', choices=['sys', 'sub'], help='which branch: system|subject')
    parser.add_argument('-s', '--switch', action='store_true', help='choose module comm_sub')
    parser.add_argument('-j', '--job', help='select one attribute')
    parser.add_argument('-k', '--subkey', help='select one key for subkeys')
    args = parser.parse_args()

    regex = re.compile('__')    # only detect it starts with lower case
    #if args.switch==False and args.mod == 'sys':
    if args.switch==False:
        mod_name = 'comment_sys'
    else:
        mod_name = 'comment_subj'
    my_module = importlib.import_module(mod_name)

    if not args.job:
        my_module.print_obj()
        if mod_name == 'comment_sys':
            print(f"\t    -s for other attributes in module 'comment_subj.py' ")
    else:
        jobs(my_module, args.job, args.subkey)

if __name__ == "__main__":
    main()
