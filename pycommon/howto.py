#!/home/joonho/anaconda3/bin/python

import argparse
import re
import importlib

### not used
def format_text(text, mod):
    if isinstance(text, str) and '{POSCAR}' in text:
        return text.format(POSCAR=getattr(mod, 'POSCAR', '{POSCAR}'))
    return text


def jobs(mod_comm, att, subkey, poscar=None):
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
    for i, key in enumerate(sel_subkeys):
        print(f"subkey-{i} {key}")          #1st att(server), 2nd sge 
        ### if second key is dict
        if isinstance(mod_comm.__dict__[att].__dict__[key], dict):   # value of 'sge' is dict
        #print("here is True")
            for k in mod_comm.__dict__[att].__dict__[key].__dict__.keys():
                #print(f" k is {k}")
                print(mod_comm.__dict__[att].__dict__[key].__dict__[k].format(POSCAR=poscar))
            #sub_keys.append(key)
        ### only 1st key (hfse2) is dict, poscar is not dict
        else:
            print(mod_comm.__dict__[att].__dict__[key].format(POSCAR=poscar))
    if not subkey:
        print(f"Use -k for subkey in {subkeys}")

    return 0

def main():
    parser = argparse.ArgumentParser(description="shows dictionary for all: work, system, package  ")
    #parser.add_argument('-m', '--mod', default='sys', choices=['sys', 'sub'], help='which branch: system|subject')
    parser.add_argument('-s', '--switch', action='store_true', help='choose module comm_sub')
    parser.add_argument('-j', '--job', help='select one attribute')
    parser.add_argument('-p', '--poscar', help='args such as POSCAR name')
    parser.add_argument('-k', '--subkey', help='select one key for subkeys')
    args = parser.parse_args()

    regex = re.compile('__')    # only detect it starts with lower case
    #if args.switch==False and args.mod == 'sys':
    if args.switch==False:
        mod_name = 'comment_subj'
    else:
        mod_name = 'comment_sys'
    my_module = importlib.import_module(mod_name)
    my_module.POSCAR = args.poscar        # module level variable: add new attribute 'POSCAR' to mod my_module
    ### same as
    #my_module.__dict__['POSCAR'] = "anything"
    ### try to pass args 
    if not args.job:
        print(my_module.__file__)
        my_module.print_obj( job = args.job, poscar=args.poscar )
        if mod_name == 'comment_subj':
            print(f"\t    -s for other attributes in module 'comment_sys.py' ")
    else:
        jobs(my_module, args.job, args.subkey, args.poscar)

if __name__ == "__main__":
    main()
