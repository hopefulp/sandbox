#!/home/joonho/anaconda3/bin/python

import argparse
import re
import comment

def jobs(att, lkeys, Lall):
    keys_ = comment.__dict__[att].__dict__.keys()
    print('---------------------------------------------')
    print(f"{att} keys:: {keys_}")
    
    if Lall:
        for key in keys_:
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
    parser.add_argument('-k', '--key_', nargs='*', help='select keys')
    parser.add_argument('-a', '--all', action='store_false', help='display all keywords for the attribute')
    args = parser.parse_args()

    regex = re.compile('[_A-Z]')
    att = [ x for x in dir(comment) if not regex.match(x) ]
    print("===ALL ATTributes in comment.py ===")
    print(f"  {att}")

    if not args.job:
        print(f"Use -j for attribute: {att}")
    else:
        jobs(args.job, args.key_, args.all)

if __name__ == "__main__":
    main()
