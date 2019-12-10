#!/home/joonho/anaconda3/bin/python

import argparse
import re
import comment

def jobs(att, lkeys,next_key, Lall):
    keyss = comment.__dict__[att].__dict__.keys()
    print('---------------------------------------------')
    print(f"{att} keys:: {keyss}")
    
    if Lall:    # Now this is running
        ### use select next key here
        for key in keyss:
            #print(key)          #1st att(server), 2nd sge 
            if isinstance(comment.__dict__[att].__dict__[key], dict):   # value of 'sge' is dict
                #print("here is True")
                for k in comment.__dict__[att].__dict__[key].__dict__.keys():
                    #print(f" k is {k}")
                    print(comment.__dict__[att].__dict__[key].__dict__[k])
            else:
                print(comment.__dict__[att].__dict__[key])
    elif lkeys:
        for key in lkeys:
            print(comment.__dict__[att].__dict__[key])
    else:
        print("Use [-a] for all keys || [-k] for some keys")

    return 0

def main():
    parser = argparse.ArgumentParser(description="shows dictionary for all: work, system, package  ")
    parser.add_argument('-j', '--job', help='select one attribute')
    parser.add_argument('-k', '--keyss', nargs='*', help='select keys')
    parser.add_argument('-k2', '--subkey', help='select keys')
    parser.add_argument('-a', '--all', action='store_false', help='display all keywords for the attribute')
    args = parser.parse_args()

    regex = re.compile('[a-z]')    # only detect it starts with lower case
    att = [ x for x in dir(comment) if regex.match(x) ]
    print("===ALL ATTributes in comment.py ===")
    print(f"  {att}")

    if not args.job:
        print(f"Use -j for attribute: {att}")
    else:
        jobs(args.job, args.keyss, args.subkey, args.all)

if __name__ == "__main__":
    main()
