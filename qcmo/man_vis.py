#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_all, MyClass, dir_classify

iqmol = MyClass()
iqmol.xyz="x, y, z-axes:: pink, green, violet-arrows"


jmol= MyClass()
jmol.orbital="orbital can be drawn separately or combined"

classobj_dict={'IQMOL':iqmol, "JMOL":jmol}

def jobs(job):
    
    if job == None:
        print(f"choose:: {classobj_dict.keys()}")
        print("#Comment: -j visual-program")
    elif job == 'IQMOL':
        attrs = vars(classobj_dict[job])
        print(attrs)
        #print(f"{item}" for item in arrts.items())
                    
    

def main():
    parser = argparse.ArgumentParser(description="display Usage for visual program ")
    parser.add_argument('-j','--job', help="select visual program such as jmol, iqmol, etc ")
    args = parser.parse_args()

    jobs(args.job)

if __name__ == "__main__":
    main()
