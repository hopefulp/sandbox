#!/home/noische/python

import sys, re, string, getopt, optparse, math, time
import bgf
import bgftools
import numpy as np
import operator
import copy
import nutils as nu
from math import *
import tqdm
import mathtools

option = ""; args = ""; bgf_file = ""; 
usage = """
PEI_findSmallestBox.py -b bgf_file -f ff_file(s) (-o out_file)
"""
version = "131216"

"""
20131216 Find the least size of the box by checking while rotation
"""
def rotate(tx, ty, tz):
    """
    Returns the rotation matrix.
    Note that tx, ty, tz are radians.
    """

    Rx = np.array([[1,0,0], [0, cos(tx), -sin(tx)], [0, sin(tx), cos(tx)]])
    Ry = np.array([[cos(ty), 0, -sin(ty)], [0, 1, 0], [sin(ty), 0, cos(ty)]])
    Rz = np.array([[cos(tz), -sin(tz), 0], [sin(tz), cos(tz), 0], [0,0,1]])

    return np.dot(Rx, np.dot(Ry, Rz))


def findSmallestBoxsize(bgf_file, out_file, ff_file, silent=False):

    # const
    PI = math.pi

    ## variables
    volume = 0.0; min_volume = 1000000000    # volumes
    angle = [];
    a = 0.0; b = 0.0; c = 0.0;    # radians
    step = 0;
    output = "";
    t1 = t2 = 0; elapsed_time = 0;


    ## Read BGF
    myBGF = bgf.BgfFile(bgf_file)

    x = np.arange(0, 2*PI, PI/9)    # x (phi)
    y = np.arange(0, 2*PI, PI/9)    # y (theta)
    z = np.arange(0, 2*PI, PI/9)    # z (rho)

    n_step = len(x) * len(y) * len(z)
    

    ## get COM
    xc, yc, zc = bgftools.getCom(myBGF, ff_file)
    for atom in myBGF.a:
        atom.x -= xc
        atom.y -= yc
        atom.z -= zc

    ## Scan
    for p in tqdm.tqdm(x, desc="Scanning", ncols=120):
        for t in y:
            for r in z:
                # rotation matrix
                U = mathtools.rotate_matrix(p, t, r)

                # rotate all atoms
                for atom in myBGF.a:
                    v = np.matrix([atom.x, atom.y, atom.z]).T
                    Uv = U * v
                    atom.x = float(Uv[0])
                    atom.y = float(Uv[1])
                    atom.z = float(Uv[2])

                #calculate the box size: check the atom position (is not related to CRYSTX)
                xlo = ylo = zlo = 1000; xhi = yhi = zhi = -1000;

                for atom in myBGF.a:
                    if atom.x < xlo:
                        xlo = atom.x
                    if atom.x > xhi:
                        xhi = atom.x
                    if atom.y < ylo:
                        ylo = atom.y
                    if atom.y > yhi:
                        yhi = atom.y
                    if atom.z < zlo:
                        zlo = atom.z
                    if atom.z > zhi:
                        zhi = atom.z

                volume = (xhi-xlo) * (yhi-ylo) * (zhi-zlo)
                
                #if its volume is the minimum so far
                if volume < min_volume:
                    # update the value
                    min_volume = volume
                    angle = [p, t, r]

                output += str(p) + "\t" + str(t) + "\t" + str(r) + "\t" + str(volume) + "\n"


    #get the rotated image from the original
    myBGF2 = bgf.BgfFile(bgf_file)
    U = rotate(angle[0], angle[1], angle[2])
    for atom in myBGF2.a:
        v = np.matrix([atom.x, atom.y, atom.z]).T
        Uv = U * v
        atom.x = float(Uv[0])
        atom.y = float(Uv[1])
        atom.z = float(Uv[2])
    
    # write output
    myBGF2.saveBGF(out_file)
    if not silent:
        print(output)
        print("Found the minimum box size: %8.3f %8.3f %8.3f %12.3f" % (angle[0], angle[1], angle[2], min_volume))

    return 1

    ### end of function


if __name__ == "__main__":
    bgf_file = ""; out_file = ""; ff_file = ""; silent = True;
    options, args = getopt.getopt(sys.argv[1:], 'hb:o:f:v', ['help', 'bgf=', 'out=', 'ff=', 'verbose'])

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
        elif option in ('-o', '--out'):
            out_file = value
        elif option in ('-f', '--ff'):
            ff_file = value
        elif option in ('-v', '--verbose'):
            silent = False
        elif option == NULL:
            print(usage)
            sys.exit(0)
    
    if out_file == "":
        out_file = bgf_file.split(".bgf")[0] + ".minsize.bgf"

    # main call
    findSmallestBoxsize(bgf_file, out_file, ff_file, silent=silent)

