#!/home/joonho/anaconda3/bin/python

import argparse
import re
import importlib

def jobs(comment, att, subkey):
    subkeys = comment.__dict__[att].__dict__.keys()
    print('---------------------------------------------')
    print(f"{att} sub_keys:: {subkeys}")
    Tsubkey = 'NO'
    sub_keys=[]
    #if isinstance(comment.__dict__[att].__dict__[keyss[0]],dict):
    ### use select next key here
    if subkey:
        subkeys=[]
        subkeys.append(subkey)
    for key in subkeys:
        #print(key)          #1st att(server), 2nd sge 
        if isinstance(comment.__dict__[att].__dict__[key], dict):   # value of 'sge' is dict
        #print("here is True")
            for k in comment.__dict__[att].__dict__[key].__dict__.keys():
                #print(f" k is {k}")
                print(comment.__dict__[att].__dict__[key].__dict__[k])
            sub_keys.append(key)
        else:
            print(comment.__dict__[att].__dict__[key])
    if not subkey:
        print(f"Use -k for sub_key in {sub_keys}")

    return 0

def main():
    parser = argparse.ArgumentParser(description="shows dictionary for all: work, system, package  ")
    parser.add_argument('-m', '--mod', default='sys', choices=['sys', 'sub'], help='which branch: system|subject')
    parser.add_argument('-s', '--switch', action='store_true', help='choose module comm_sub')
    parser.add_argument('-j', '--job', help='select one attribute')
    parser.add_argument('-k', '--subkey', help='select one key for subkeys')
    args = parser.parse_args()

    regex = re.compile('[a-z]')    # only detect it starts with lower case
    if args.switch==False and args.mod == 'sys':
        mod_name = 'comm_sys'
    else:
        mod_name = 'comm_sub'
    my_module = importlib.import_module(mod_name)
    att = [ x for x in dir(my_module) if regex.match(x) ]
    print(f"===ALL ATTributes in \"{mod_name}.py\" ===")
    print(f"  {att}")

    if not args.job:
        #print(f"Use -j for attribute: {att}")
        print(f"Use -j for attribute")
        if mod_name == 'comm_sys':
            print(f"Use -m sub or -s for 'comm_sub.py' if attribute you find is not here")
    else:
        jobs(my_module, args.job, args.subkey)

if __name__ == "__main__":
    main()
