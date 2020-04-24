#!/home/joonho/anaconda3/bin/python

"""
The script converts XYZ file to lammps data in atomic format for Tersoff potential.
i.e. type(element) x y z -> id type(number) x y z
"""

import csv
import sys
import os
import warnings
import pandas as pd
import argparse
from common import f_root, f_ext

atom_mass = {'H': 1.01, 'He': 4.00, 'Li': 6.94, 'Be': 9.01, 'B': 10.81, 'C': 12.01, 'N': 14.01, 'O': 16.00, 'F': 19.00,
             'Ne': 20.18, 'Na': 22.99, 'Mg': 24.31, 'Al': 26.98, 'Si': 28.09, 'P': 30.97, 'S': 32.07, 'Cl': 35.45,
             'Ar': 39.95, 'K': 39.10, 'Ca': 40.08, 'Sc': 44.96, 'Ti': 47.87, 'V': 50.94, 'Cr': 52.00, 'Mn': 54.94,
             'Fe': 55.85, 'Co': 58.93, 'Ni': 58.69, 'Cu': 63.55, 'Zn': 65.39, 'Ga': 69.72, 'Ge': 72.61, 'As': 74.92,
             'Se': 78.96, 'Br': 79.90, 'Kr': 83.80, 'Rb': 85.47, 'Sr': 87.62, 'Y': 88.91, 'Zr': 91.22, 'Nb': 92.91,
             'Mo': 95.94, 'Tc': 98.00, 'Ru': 101.07, 'Rh': 102.91, 'Pd': 106.42, 'Ag': 107.87, 'Cd': 112.41,
             'In': 114.82, 'Sn': 118.71, 'Sb': 121.76, 'Te': 127.60, 'I': 126.90, 'Xe': 131.29, 'Cs': 132.91,
             'Ba': 137.33, 'La': 138.91, 'Ce': 140.12, 'Pr': 140.91, 'Nd': 144.24, 'Pm': 145.00, 'Sm': 150.36,
             'Eu': 151.96, 'Gd': 157.25, 'Tb': 158.93, 'Dy': 162.50, 'Ho': 164.93, 'Er': 167.26, 'Tm': 168.93,
             'Yb': 173.04, 'Lu': 174.97, 'Hf': 178.49, 'Ta': 180.95, 'W': 183.84, 'Re': 186.21, 'Os': 190.23,
             'Ir': 192.22, 'Pt': 195.08, 'Au': 196.97, 'Hg': 200.59, 'Tl': 204.38, 'Pb': 207.2, 'Bi': 208.98,
             'Po': 209.00, 'At': 210.00, 'Rn': 222.00, 'Fr': 223.00, 'Ra': 226.00, 'Ac': 227.00, 'Th': 232.04,
             'Pa': 231.04, 'U': 238.03, 'Np': 237.00, 'Pu': 244.00, 'Am': 243.00, 'Cm': 247.00, 'Bk': 247.00,
             'Cf': 251.00, 'Es': 252.00, 'Fm': 257.00, 'Md': 258.00, 'No': 259.00, 'Lr': 262.00, 'Rf': 261.00,
             'Db': 262.00, 'Sg': 266.00, 'Bh': 264.00, 'Hs': 269.00, 'Mt': 268.00}

def xyz2lmp(xyz_file, out_file, lattice_a):

    yes_bgf = False

    # open xyz file
    print("Reading the file..")
    molecule = pd.read_table(xyz_file, skiprows=2, delim_whitespace=True, names=['atom', 'x', 'y', 'z'])

    # check whether similar BGF file exists
    bgf_file = xyz_file.split(".xyz")[0] + ".bgf"
    #if os.path.exists(bgf_file):
    #    print("A BGF file similar to the xyz file found (%s). Will use CRYSTX info from BGF." % bgf_file)
    #    yes_bgf = True
    #    from bgf.BGF import *
    #    bgf_file = BgfFile(bgf_file)

    # extract atom types
    print("Processing..")
    atoms = {}
    for index, i in enumerate(sorted(list(set(list(molecule['atom']))))):
        atoms[i] = index + 1

    # create an atom_id column
    # molecule['atom'] += 1
    # molecule = molecule[['atom', 'x', 'y', 'z']]

    atom_id = [atoms[i] for i in molecule['atom']]
    atom_id_series = pd.Series(atom_id, name='atom_id')
    # molecule.insert(loc=1, column='atom_id', value=atom_id_series)
    molecule.insert(loc=1, column='atom_id', value=atom_id)
    del molecule['atom']
    # molecule.insert(0, 'atom_id', atom_id_series)
    molecule.index += 1

    # atom numbers and types
    header = "# LAMMPS Tersoff Potential Input converted from %s\n\n" % xyz_file
    header += "%11s atoms\n%11s atom types\n\n" % (len(molecule), len(atoms))

    # simulation box
    if yes_bgf:
        header += "%16.8f %16.8f    xlo xhi\n" % (0, bgf_file.CRYSTX[0])
        header += "%16.8f %16.8f    ylo yhi\n" % (0, bgf_file.CRYSTX[1])
        header += "%16.8f %16.8f    zlo zhi\n" % (0, bgf_file.CRYSTX[2])
    else:
        header += "%16.8f %16.8f    xlo xhi\n" % (0, lattice_a)
        header += "%16.8f %16.8f    ylo yhi\n" % (0, lattice_a)
        header += "%16.8f %16.8f    zlo zhi\n" % (0, lattice_a)
        #print("WARNING: PLEASE DO NOT FORGET TO FILL OUT XLO XHI YLO YHI ZLO ZHI TERMS.")

    # atom masses
    header += "\nMasses\n\n"
    for key, value in atoms.items():
        header += "%s    %.3f  # %s\n" % (value, atom_mass[key], key)

    # atom coordinates
    header += "\nAtoms\n\n"
    with open(out_file, "wt") as f:
        f.write(header)

    # save
    molecule.to_csv(path_or_buf=out_file, sep="\t", float_format="%11.5f", index=True, header=False, mode="a")

    # print("Postprocessing..")
    # with open(out_file, 'rt') as f:
    #     line = f.readlines()
    #
    # with open(out_file, 'wt') as f:
    #     for l in line:
    #         f.write(l.replace("\t", "    "))
    # molecule.to_string(buf=out_file, float_format="%11.5f", index_names=True, header=False)

    print("Done. File save to %s !" % out_file)
    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="from xyz file to make lammps input data file")
    parser.add_argument('infile',  help="xyz input file")
    parser.add_argument('-o', '--outfile', help="lammps input data file")
    parser.add_argument('-a', '--lattice_a', default=10.0, type=float, help="lattice vector of a")
    args = parser.parse_args()

    if f_ext(args.infile) != "xyz":
        print("input xyz file")
        sys.exit(1)

    if not args.outfile:
        outfile = "data."+f_root(args.infile)
    else:
        outfile = args.outfile

    xyz2lmp(args.infile, outfile, args.lattice_a)
