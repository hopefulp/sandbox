#!/home/noische/program/python27/bin/python

# python modules
import sys
import os
import copy
import re
from types import *

import networkx as nx

version = '110819'

masses = {}

def get_mass(dat_file):

    parse = [line.rstrip() for line in open(dat_file, 'r')]

    for i in parse:
        if "atom types" in i:
            n_atom_types = int(re.split('\s*', i.lstrip())[0])

    start = parse.index("Masses") + 2
    slice = parse[start:start + n_atom_types]

    for i in slice:
        line = re.split('\s*', i.strip())
        masses[int(line[0])] = float(line[1])

    return masses


def get_connection(dat_file):

    parse = [line.rstrip() for line in open(dat_file, 'r')]

    for i in parse:
        if "bonds" in i:
            n_bonds = int(re.split('\s*', i.lstrip())[0])

    start = parse.index("Bonds") + 2
    slice = parse[start:start + n_bonds]

    bonds = []
    for i in slice:
        line = i.split()
        atom1 = int(line[2]); atom2 = int(line[3])
        bonds.append([atom1, atom2])

    g = nx.Graph(bonds)
    connection = sorted([list(i) for i in nx.connected_components(g)], key=lambda x:x[0])

    return connection


def get_atom_types(dat_file):

    parse = [line.rstrip() for line in open(dat_file, 'r')]

    for i in parse:
        if "atoms" in i:
            n_atoms = int(re.split('\s*', i.lstrip())[0])

    start = parse.index("Atoms") + 2
    slice = parse[start:start + n_atoms]

    atoms = {}
    for i in slice:
        line = i.lstrip().split()
        atom_no = int(line[0]); atom_type = int(line[2])
        atoms[atom_no] = atom_type

    return atoms

