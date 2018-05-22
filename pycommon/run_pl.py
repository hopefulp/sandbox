#!/usr/bin/python

import argparse
import os
from common import *


def run_perl(pl, T_matching, pattern, exe):
    pwd = os.getcwd()

    print pattern
    l_file = get_files_pattern(T_matching, pattern, pwd)
    
    for file in l_file:
        sh_comm = pl + file
        print sh_comm
        if exe:
            os.system(sh_comm)

    return (0)        

def main():
    parser = argparse.ArgumentParser(description="run perl script to many files")
    parser.add_argument("plscript", help="perl script")
    parser.add_argument('-t', '--style', choices=['p', 's', 'm'], help="matching method, p[refix], s[uffix], m for searching")
    parser.add_argument('-e', '--exe', action="store_true", help="execute perl command")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-p', '--prefix', nargs='+', help="prefix exclusive to suffix and search")
    group.add_argument('-s', '--suffix', nargs='+', help="suffix exclusive to prefix and search")
    group.add_argument('-m', '--matching', nargs='+', help="matching(re.search) exclusive to prefix and suffix")
    args = parser.parse_args()

    if args.prefix:
        matching=args.prefix
        m_tag = 'p'
    elif args.suffix:
        matching=args.suffix
        m_tag = 's'
    elif args.matching:
        matching=args.matching
        m_tag = 'm'
    else:
        print "matching should be given"
        return 1

    run_perl(args.plscript, m_tag, matching, args.exe)
    return 0

if __name__ == "__main__":
    main()
