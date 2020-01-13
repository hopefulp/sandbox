#!/home/joonho/anaconda3/bin/python

import argparse
import re
import comment
import comment_subj

def jobs(att, subkey):
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
    parser.add_argument('-c', '--classify', default='general', choices=['general', 'subjects'], help='which branch: general|subject')
    parser.add_argument('-j', '--job', help='select one attribute')
    parser.add_argument('-k', '--subkey', help='select one key for subkeys')
    args = parser.parse_args()

    regex = re.compile('[a-z]')    # only detect it starts with lower case
    att = [ x for x in dir(comment) if regex.match(x) ]
    print("===ALL ATTributes in comment.py ===")
    print(f"  {att}")

    if not args.job:
        print(f"Use -j for attribute: {att}")
    else:
        jobs(args.job, args.subkey)

if __name__ == "__main__":
    main()
