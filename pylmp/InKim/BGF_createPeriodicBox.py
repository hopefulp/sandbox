#!/home/noische/python

import sys
import re
import string
import getopt
import optparse
import math
import time

from os import popen
import bgf

#-----------------
# update the coordinate in the original BGF file from LAMMPS trajectory file
#_________________
def createPeriodicBox(bgf_file, out_file, str_pbc, silent=True):
    boxsize = [0, 0, 0, 0, 0, 0]

    # open bgf
    if isinstance(bgf_file, bgf.BgfFile):
        myBGF = bgf_file
    else:
        if not silent: print("reading " + bgf_file + " ..")
        myBGF = bgf.BgfFile(bgf_file)

    myBGF.PERIOD = "111"
    myBGF.AXES = "ZYX"
    myBGF.SGNAME = "P 1                  1    1"
    myBGF.CELLS = [-1, 1, -1, 1, -1, 1]

    if str_pbc == "":
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
            
        ### Size of the box
        a = []
        a.append(xhi-xlo+0)
        a.append(yhi-ylo+0)
        a.append(zhi-zlo+0)
        a.append(90.0)
        a.append(90.0)
        a.append(90.0)
        myBGF.CRYSTX = a
    else:
        pbc = str_pbc.split()
        try:
            pbc = [float(i) for i in pbc]
        except:
            nu.die("Wrong PBC conditions given.")

        if len(pbc) == 3:
            myBGF.CRYSTX = pbc + [90.0, 90.0, 90.0]
        elif len(pbc) == 6:
            myBGF.CRYSTX = pbc
        else:
            nu.die("Wrong PBC conditions given.")

    # save
    if isinstance(out_file, str):
        if not silent: print("saving information to " + out_file + " ..")
        myBGF.saveBGF(out_file)
        return 1;
    else:
        return myBGF;


if __name__ == "__main__":

    option = ""; args = ""; bgf_file = ""; out_file = ""; pbc = "";
    usage = """
Usage: createPeriodicBox.py -b bgf_file -o out_file
    Creates a periodic box to a BGF file.
    It is highly recommended to center your BGF with ~tpascal/scripts/centerbgf.pl
    """

    if len(sys.argv) < 2:
        print(usage)
        sys.exit(0)

    print(sys.argv)

    options, args = getopt.getopt(sys.argv[1:], 'hb:o:p:', ['help','bgf=','out=', 'pbc='])
    for option, value in options:
        if option in ('-h', '--help'):
            print(usage); sys.exit(0)
        elif option in ('-b', '--bgf'):
            bgf_file = value
        elif option in ('-o', '--out'):
            out_file = value
        elif option in ('-p', '--pbs'):
            pbc = str(value)
        elif option in NULL:
             print(usage); sys.exit(0)

    # main call
    createPeriodicBox(bgf_file, out_file, pbc)
