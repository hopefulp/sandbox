#!/home/noische/program/python27/bin/python

import sys
import os
import string
import getopt
import pprint
import random
import cutils as cu
import nutils as nu
import bgftools

version = '110101'

#-------------------------------------
#
# - Main Function
#
"""if __name__ == '__main__':

	bgftools.renumberMolecules("methanol.box.dreiding.bgf", "methanol.box.dreiding.mod.bgf", False)
"""

def test(*args, **kwargs):
    print(args)
    print(kwargs)
    if not 'mark' in kwargs:
        mark = True

    print(mark)

print("test 1")
test()
print('test 2')
test(mark='a')
