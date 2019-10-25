#!/home/noische/Enthought/Canopy_64bit/User/bin/python
### PRACTICE ###

import sys

def read_big_file(file):

    def get_line(fname):
        with open(fname) as f:
            for line in f:
                yield line

    dump = get_line(file)
    
    for i in dump:
        print i.strip('\n')

read_big_file(sys.argv[1])
