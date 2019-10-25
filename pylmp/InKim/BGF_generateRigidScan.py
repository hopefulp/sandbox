#!/home/noische/python

import sys
import getopt
import copy

import numpy as np
import tqdm

import bgf
import bgftools
import nutils as nu

usage = """
Generates BGF file with changed distances between two molecules.
This is useful to generate input structures for Rigid Coordinate Scan.

Usage: %s -b bgf_file -1 atom1 -2 atom2 -o out_file -d "distances"
    -b bgf_file: self-descriptive.
    -1 atom1: atom number to fix.
    -2 atom2: atom number to change distance. should be in another atom to atom1.
    -o out_file: self-descriptive.
    -d distances: self-descriptive. optional.
                "-1.0 1.0 2.0" will generate three structures with those distances changed.
                defaults: "-1.0 0.0 1.0 2.0 ... 10.0"
    
Please report any bugs to in.kim@kaist.ac.kr

""" % sys.argv[0]

def generate_rcs(bgf_file, atom1, atom2, out_file, dist=[], silent=False):
    # init
    if not dist:
        dist = np.arange(-1.0, 10.0, 1.0)

    # open
    if isinstance(bgf_file, bgf.BgfFile):
        myBGF = bgf_file
    else:
        if not silent: print("reading " + bgf_file + " ..")
        myBGF = bgf.BgfFile(bgf_file)

    a1 = myBGF.a[myBGF.a2i[atom1]];
    a2 = myBGF.a[myBGF.a2i[atom2]];

    # check whether a1 and a2 are in a same molecule // the script cannot continue if yes
    mol1 = []; mol2 = []
    dummy = bgftools.getmolecule(myBGF, a1, mol1)
    dummy = bgftools.getmolecule(myBGF, a2, mol2)
    if mol1 == mol2:
        nu.die("atom1 and atom2 are in a same molecule: cannot perform a rigid body scan within a same molecule.")

    mybgf2 = copy.deepcopy(myBGF)

    # translate the coordinate: atom1 to origin
    v1 = (a1.x, a1.y, a1.z) # position of atom1
    for atom in mybgf2.a:
        atom.x -= v1[0]
        atom.y -= v1[1]
        atom.z -= v1[2]

    # find the rotation matrix: v21 to x axis
    v21 = [ (a2.x - a1.x), (a2.y - a1.y), (a2.z - a1.z) ]
    u1 = v21 / np.linalg.norm(v21)    ##
    v2 = [u1[1], -u1[0], 0]; u2 = v2 / np.linalg.norm(v2)    ##
    v3 = np.cross(u1, u2); u3 = v3 / np.linalg.norm(v3)    ##

    U = np.array([[u1[0], u1[1], u1[2]], [u2[0], u2[1], u2[2]], [u3[0], u3[1], u3[2]]])
    Uinv = np.linalg.inv(U)

    for d in tqdm.tqdm(dist, ncols=120, desc="Creating RCS structures"):
        mybgf2 = copy.deepcopy(myBGF)

        # translate the coordinate: atom1 to origin
        v1 = (a1.x, a1.y, a1.z) # position of atom1
        for atom in mybgf2.a:
            atom.x -= v1[0]
            atom.y -= v1[1]
            atom.z -= v1[2]

        # rotate all atoms
        for atom in mybgf2.a:
            a = np.matrix([atom.x, atom.y, atom.z]).T
            b = U*a
            atom.x = float(b[0])
            atom.y = float(b[1])
            atom.z = float(b[2])

        # move mol2 by d
        for ano in mol2:
            atom = mybgf2.getAtom(ano)
            atom.x += d

        # re-rotate all atoms 
        for atom in mybgf2.a:
            a = np.matrix([atom.x, atom.y, atom.z]).T
            b = Uinv*a
            atom.x = float(b[0])
            atom.y = float(b[1])
            atom.z = float(b[2])

        ## save
        fname = out_file.split(".bgf")[0] + "." + str("{0:.2f}".format(d)).replace(".", "p") + ".bgf"
        mybgf2.saveBGF(fname)

    ### end of function


if __name__ == "__main__":
    bgf_file = ""; trj_file = ""; atom1 = 0; atom2 = 0; out_file = ""; distance = ""

    options, args = getopt.getopt(sys.argv[1:], 'hb:1:2:o:d:', ['help', 'bgf=', 'atom1=', 'atom2=', 'out=', 'distance='])

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
        elif option in ('-d', '--distance'):
            distance = str(value)
        elif option == NULL:
            print(usage)
            sys.exit(0)
    

    # scan distances:
    dist = [float(i) for i in distance.split(' ')]

    # main call
    generate_rcs(bgf_file, atom1, atom2, out_file, dist=dist)
