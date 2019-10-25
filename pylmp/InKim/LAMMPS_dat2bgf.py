#!/home/noische/python

import sys, re, string, getopt, optparse, math, time, pprint, os
from os import popen
import bgf
import bgftools
import numpy
import operator
import copy
import nutils as nu

usage = """
LAMMPS_dat2bgf.py -b bgf_file -d dat_file
"""
version = "160225"

#_________________
# 130614: applies box information from trajectory
#_________________

def dat2bgf(bgf_file, dat_file, out_file):

    # init
    mybgf = bgf.BgfFile(bgf_file)
    
    f_dat_file = open(dat_file)
    temp = f_dat_file.read().split('\n')
    f_dat_file.close()

    temp = [i.partition('#')[0].rstrip() for i in temp if i != ""]

    # PBC
    for i in temp:
        if 'xlo' in i:
            line = i.split()
            mybgf.CRYSTX[0] = float(line[1]) - float(line[0])
        if 'ylo' in i:
            line = i.split()
            mybgf.CRYSTX[1] = float(line[1]) - float(line[0])
        if 'zlo' in i:
            line = i.split()
            mybgf.CRYSTX[2] = float(line[1]) - float(line[0])

    # Coordinates
    dat_atom_start = temp.index('Atoms')
    dat_atom_end = temp.index('Bonds')
    if 'Velocities' in temp:
        dat_atom_end = temp.index('Velocities') # if data file is from write_restart command, Velocities are also written in the file.:w

    temp = temp[dat_atom_start + 1:dat_atom_end]

    n_dat_atoms = len(temp)
    n_bgf_atoms = len(mybgf.a)

    if n_dat_atoms != n_bgf_atoms:
        nu.die("Number of atoms mismatch between bgf and LAMMPS data file.")

    coords = [] # coords = id type molid charge x y z ix iy iz
    for i in temp:
        line = i.split()
        coords.append(line)

    for i in coords:
        id = int(i[0])
        atom = mybgf.a[mybgf.a2i[id]]   # getAtom
        atom.x = float(i[4])
        atom.y = float(i[5])
        atom.z = float(i[6])

    if out_file == "": out_file = bgf_file.split(".bgf")[0] + "_mod.bgf"
    mybgf.REMARK.append("Coordinates updated on %s" % time.asctime(time.gmtime()))
    mybgf.REMARK.append("Coordinates updated with data file %s" % os.path.abspath(dat_file))
    mybgf.saveBGF(out_file)


    print('Done. Check %s ' % out_file)
    return 1

    ### end of function


if __name__ == "__main__":
    bgf_file = ""; dat_file = ""; out_file = "";

    options, args = getopt.getopt(sys.argv[1:], 'hb:d:o:', ['help', 'bgf=', 'dat=', 'out='])

    if len(sys.argv) < 2:
        print(usage)
        sys.exit(0)

    print "Requested options: " + str(options)

    for option, value in options:
        if option in ('-h', '--help'):
            print(usage)
            sys.exit(0)
        elif option in ('-b', '--bgf'):
            bgf_file = value
        elif option in ('-d', '--dat'):
            dat_file = value
        elif option in ('-o', '--out'):
            out_file = value
        elif option == NULL:
            print(usage)
            sys.exit(0)
    
    # main call
    dat2bgf(bgf_file, dat_file, out_file)
