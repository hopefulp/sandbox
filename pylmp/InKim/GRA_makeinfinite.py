#!/home/noische/Enthought/Canopy_64bit/User/bin/python

"""
GRA_makeinfinite.py

Connect C atoms to make infinite graphene sheet
Most functions are from CNT.py
"""

import os, sys, copy
import numpy as np
import bgf, bgftools
import nutils as nu
from tqdm import tqdm

version = "150204"
usage = "GRA_makeinfinite.py bgf_file out_file"

def check_graphene_atom(atom):
    if not "GR" in atom.rName:
        return False
    if len(atom.CONECT) == 3:
        return False

    return True

def make_infinite(mybgf):
    def calc_bond_dist(pbc=False):
        dist = [];
        for atom in mybgf.a:
            if not "C_2G" in atom.ffType:
                continue
            for ano in atom.CONECT:
                atom2 = mybgf.getAtom(ano)
                if not "C_2G" in atom2.ffType:
                    continue

                if pbc:
                    dist.append(bgf.pbc_distance(atom, atom2, mybgf.CRYSTX[:3]))
                else:
                    dist.append(bgf.distance(atom, atom2))

        return np.average(dist), np.std(dist)

    distance, _ = calc_bond_dist()

    # TODO: need to distinguish armchair and zigzag boundaries
    minmax_x = bgftools.atoms_minmax(mybgf, "x", selection = "'C_2G' in atom.ffType")
    minmax_y = bgftools.atoms_minmax(mybgf, "y", selection = "'C_2G' in atom.ffType")

    pbc_x = (minmax_x[1] - minmax_x[0]) + distance * np.cos(np.radians(30))    # in 5nm x 5nm graphene, x boundaries are armchair
    pbc_y = (minmax_y[1] - minmax_y[0]) + distance    # in 5nm x 5nm graphene, y boundaries are zigzag

    mybgf.PERIOD = "111"
    mybgf.AXES = "ZYX"
    mybgf.SGNAME = "P 1                  1    1"
    mybgf.CELLS = [-1, 1, -1, 1, -1, 1]
    mybgf.CRYSTX = [pbc_x, pbc_y, 10.0, 90.0, 90.0, 90.0]
    print("Assigning new pbc: %s" % mybgf.CRYSTX)

    # init
    aNo_pair = [];
    pbc = mybgf.CRYSTX[:3]

    edge_atoms = []

    # loop
    for atom in tqdm(mybgf.a, ncols=80, miniters=100, desc="Analyzing atoms"):
        if not "C" in atom.ffType:
            continue
        if len(atom.CONECT) == 3:
            continue

        x = np.array([atom.x, atom.y, atom.z])
        min_atom_aNo = 100000;    # placeholder
        min_atom_d = 10000.0;

        for atom2 in mybgf.a:
            if len(atom2.CONECT) == 3:
                continue

            if atom.aNo != atom2.aNo:
                if (not atom2.aNo in atom.CONECT) and len(atom2.CONECT) != 3 and (atom.rName == atom2.rName):
                    y = np.array([atom2.x, atom2.y, atom2.z])
                    d = nu.pbc_dist(x, y, pbc)
                    if d < min_atom_d:
                        min_atom_aNo = atom2.aNo
                        min_atom_d = d
                    
        #aNo_pair.append([atom.aNo, min_atom_aNo])
        if min_atom_aNo == 100000:
            nu.die("Failed to find the closest atom for atom number " + str(atom.aNo))

        mybgf.connectAtoms(atom.aNo, min_atom_aNo)

    #for i in aNo_pair:
    #    mybgf.connectAtoms(i[0], i[1])

    check_connectivity();

    # check bond distance stdev
    _, stdev = calc_bond_dist(pbc=True)
    if stdev <= 1e5:
        print("The structure seems to have regular distance over pbc.")
    else:
        nu.die("PBC assignment does not give a regular distances over pbc.")

    if out_file == "":
        filename = bgf_file.split(".bgf")[0] + ".infinite.bgf"
    else:
        filename = out_file
    value = raw_input("Filename to save [" + filename + "]? ") or filename
    mybgf.saveBGF(value)

    ### end of make_infinite


def check_connectivity():
    # connectivity check-up
    print("Checking nanotube connectivity..")
    for atom in mybgf.a:
        if not "C" in atom.ffType:
            continue
        if len(atom.CONECT) != 3:
            nu.warn("Defect on nanotube connection found: " + str(atom.CONECTline()))

if len(sys.argv) < 2:
    print usage
    sys.exit(0)

bgf_file = sys.argv[1]
out_file = sys.argv[2]

#nu.warn(bgf_file)

#pbc = []
mybgf = bgf.BgfFile(bgf_file)

make_infinite(mybgf)
