#!/home/noische/program/python27/bin/python
"""
addIons.py
Original: Dec 28 2011 In Kim
151216: removed diameter restraint
"""

# Python Modules
import sys
import os
import string
import random
import time
import getopt

# Custom Modules
sys.path.append("/home/noische/scripts")
sys.path.append("/home/noische/script")
import bgf
import bgftools
import nutils as nu

# Globals
version = '151216'

def addIons(bgf_file, out_file, region, number, silent=False):
    """
def addIons():
    Write what this function does.

Function Parameters:
    -b bgf_file    A string of filename or BgfFile class.
    -o out_file    A string of filename or BgfFile class.
    -n number
    -r region: "xlo xhi ylo yhi zlo zhi"
    """
    # open BGF
    if isinstance(bgf_file, bgf.BgfFile):
        myBGF = bgf_file
    else:
        if not silent: print("opening bgf file.. " + str(bgf_file))
        myBGF = bgf.BgfFile(bgf_file)

    boxsize = [ myBGF.CRYSTX[0], myBGF.CRYSTX[1], myBGF.CRYSTX[2] ]

    ### find the last HETATM
    aNo_lastHETAtom = 0;
    for i in myBGF.a:
        aNo = i.aNo
        if i.aTag == 1:
            if aNo > aNo_lastHETAtom:
                aNo_lastHETAtom = aNo

    for i in range(number):
        
        """
        randx = random.uniform(0, boxsize[0])
        randy = random.uniform(0, boxsize[1])
        randz = random.uniform(0, boxsize[2])
        """
        randx = random.uniform(region[0], region[1])
        randy = random.uniform(region[2], region[3])
        randz = random.uniform(region[4], region[5])
    
        # Na+ ion
        atom_Na = bgf.BgfAtom()
        atom_Na.x = randx
        atom_Na.y = randy
        atom_Na.z = randz
        atom_Na.aTag = 1    # HETATM
        atom_Na.ffType = "Na"
        atom_Na.rName = "ION"
        atom_Na.aName = "Na"
        atom_Na.charge = 1.00
        atom_Na.rNo = 0
        atom_Na.chain = "X"
        
        # Cl- ion
        atom_Cl = bgf.BgfAtom()
        atom_Cl.x = randx + 2
        atom_Cl.y = randy + 2
        atom_Cl.z = randz + 2
        atom_Cl.aTag = 1    # HETATM
        atom_Cl.ffType = "Cl"
        atom_Cl.rName = "ION"
        atom_Cl.aName = "Cl"
        atom_Cl.charge = -1.00
        atom_Cl.rNo = 0
        atom_Cl.chain = "X"

        myBGF.addAtom(atom_Na, myBGF.a2i[aNo_lastHETAtom + 1])
        myBGF.addAtom(atom_Cl, myBGF.a2i[aNo_lastHETAtom + 2])

    myBGF.renumber()

    # save BGF
    if isinstance(out_file, str):
        if not silent: print("Saving information to " + out_file + " ..")
        myBGF.saveBGF(out_file)
        return 1;
    else:
        return myBGF;


    ### end of addIons


if __name__ == '__main__':

    option = ""; args = ""; bgf_file = ""; size = 0.0; out_file = ""; number = 0;
    region = []
    usage = """
Usage: addIons.py -b "bgfFiles" -s size -n monomers -o output

Options are:
    -b    Input BGF file.
    -o    Output BGF file.
    -n    Number of Na+ and Cl-
    """

    if len(sys.argv) < 2:
        print(usage); sys.exit(0)

    options, args = getopt.getopt(sys.argv[1:], 'hb:o:n:r:', ['help','bgf=','output=','number=', 'region='])
    for option, value in options:
        if option in ('-h', '--help'):
            print usage; sys.exit(0)
        elif option in ('-b', '--bgf'):
            bgf_file = value
        elif option in ('-o', '--output'):
            out_file = value
        elif option in ('-n', '--number'):
            number = int(value)
        elif option in ('-r', '--region'):
            region = str.split(value)
        elif option in (''):
            print(usage); sys.exit(0)

    # default settings
    if not out_file: out_file = os.path.basename(bgf_file).split(".bgf")[0] + "_ion" + ".bgf"
    region = [float(x) for x in region]

    addIons(bgf_file, out_file, region, number, silent=False)

