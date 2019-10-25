#!/opt/applic/epd/bin/python

import sys, re, string, getopt, optparse, math, os

import numpy as np

import bgf
import bgftools
import nutils as nu

usage = """
alignBgf.py -b bgf_file -1 atom1 -2 atom2 -o out_file
"""

def alignBGF(bgf_file, atom1, atom2, out_file, silent=False):

    # open
    if isinstance(bgf_file, bgf.BgfFile):
        myBGF = bgf_file
    else:
        if not silent: print("reading " + bgf_file + " ..")
        myBGF = bgf.BgfFile(bgf_file)

    a1 = myBGF.a[myBGF.a2i[atom1]];
    a2 = myBGF.a[myBGF.a2i[atom2]];

    v1 = (a1.x, a1.y, a1.z)

    # move atom1 to origin
    for atom in myBGF.a:
        atom.x -= v1[0]
        atom.y -= v1[1]
        atom.z -= v1[2]

    v21 = [ (a2.x - a1.x), (a2.y - a1.y), (a2.z - a1.z) ]

    # rotate v21 to x axis
    u1 = v21 / np.linalg.norm(v21)    ##

    v2 = [u1[1], -u1[0], 0]
    u2 = v2 / np.linalg.norm(v2)    ##

    v3 = np.cross(u1, u2)
    u3 = v3 / np.linalg.norm(v3)    ##

    U = np.array([[u1[0], u1[1], u1[2]], [u2[0], u2[1], u2[2]], [u3[0], u3[1], u3[2]]])

    # rotate all atoms
    for atom in myBGF.a:
        a = np.matrix([atom.x, atom.y, atom.z]).T
        b = U*a
        atom.x = float(b[0])
        atom.y = float(b[1])
        atom.z = float(b[2])

    ## save
    if isinstance(out_file, str):
        if not silent: print("saving information to " + out_file + " ..")
        myBGF.saveBGF(out_file)
        return 1;
    else:
        return myBGF;

    ### end of function


if __name__ == "__main__":
    bgf_file = ""; trj_file = ""; atom1 = 0; atom2 = 0; out_file = ""

    options, args = getopt.getopt(sys.argv[1:], 'hb:1:2:o:', ['help', 'bgf=', 'atom1=', 'atom2=', 'out='])

    if len(sys.argv) < 2:
        print(usage)
        sys.exit(0)

    print("Requested options: " + str(options))

    for option, value in options:
        if option in ('-h', '--help'):
            print(usage)
            sys.exit(0)
        elif option in ('-b', '--bgf'):
            bgf_file = value
        elif option in ('-1', '--atom1'):
            atom1 = int(value)
        elif option in ('-2', '--atom2'):
            atom2 = int(value)
        elif option in ('-o', '--out'):
            out_file = value
        elif option == NULL:
            print(usage)
            sys.exit(0)
    
    # main call
    alignBGF(bgf_file, atom1, atom2, out_file)
