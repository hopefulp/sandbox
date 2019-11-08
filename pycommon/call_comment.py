#!/home/joonho/anaconda3/bin/python

import argparse
import re
import comment

def jobs(att, lkeys):
    keys_ = comment.__dict__[att].__dict__.keys()
    print('---------------------------------------------')
    print(f"{att} keys:: {keys_}")
    
    if lkeys:
        if lkeys[0] == 'a':
            for key in keys_:
                print(comment.__dict__[att].__dict__[key])
        else:
            for key in lkeys:
                print(comment.__dict__[att].__dict__[key])
    else:
        print("Use -k a to display all keys")

    return 0

def main():
    parser = argparse.ArgumentParser(description="shows dictionary for all: work, system, package  ")
    parser.add_argument('-a', '--att', help='select one attribute')
    parser.add_argument('-k', '--key_', nargs='*', help='select keys')
    args = parser.parse_args()

    regex = re.compile('[_A-Z]')
    att = [ x for x in dir(comment) if not regex.match(x) ]
    print("===ALL ATTributes in comment.py ===")
    print(f"  {att}")

    if not args.att:
        print("Use -a for attribute")
    else:
        jobs(args.att, args.key_)

if __name__ == "__main__":
    main()
