#!/usr/bin/python

import argparse
import re
import os
import sys
from common import *

def cmd(d_ref, d_target, rm, d_trash):

    ld_target   = os.listdir(d_target)
    ld_ref      = os.listdir(d_ref)
    #print ld_target
    #print ld_ref
    for f in ld_target:
        fname = d_target + '/' + f
        if os.path.isfile(fname):
            #print f
            if f in ld_ref:
                comm = "cmp %s/%s %s/%s" % (d_target, f, d_ref, f)
                answer = os.system(comm) # if same, it returns 0
                if not answer:
                    #print "%s will be deleted in %s" % (f, d_target)
                    if os.path.isdir(d_trash):
                        comm = "mv %s/%s %s" % (d_target, f, d_trash)
                        print comm
                        if rm:
                            os.system(comm)
    return 0 


def main():
    parser = argparse.ArgumentParser(description='directory comparison ')
    parser.add_argument( 'dir_ref', help='reference directory')
    parser.add_argument( 'dir_target', help='target directory')
    parser.add_argument( '-r', '--remove', action='store_true', help='action remove')
    parser.add_argument( '-t', '--trash_dir', default="/scratch/Trash", help='trash box')
    args = parser.parse_args()

    cmd(args.dir_ref, args.dir_target, args.remove, args.trash_dir)
    return 0

if __name__ == "__main__":
    main()
