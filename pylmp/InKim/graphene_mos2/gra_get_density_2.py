#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys
import getopt
import optparse
import cPickle as pickle

import numpy as np
from numpy import arccos
from numpy.linalg import norm
import tqdm

import bgf
import bgftools as bt
import lammpstrj as lt
#from lammpstrj import *
import nutils as nu

from lammps_analyze_trj_2 import *

def count_layer(bgf_file, trj_file, ff_file='', out_file=''):
    '''counts number of O atoms in each layer.
    '''
    # variables
    result = dict()

    # 1. Load BGF
    mybgf = bgf.BgfFile(bgf_file)
    N_BGF_ATOMS = len(mybgf.a)
    r_vdw_C = 3.38383824/2

    # 2. Read LAMMPS Trajectory
    mytrj = lt.lammpstrj(trj_file)
    timesteps = mytrj.load()
    N_ATOMS = mytrj.natoms[0]
    if N_BGF_ATOMS != N_ATOMS:
        nu.die("Number of atoms in trajectory file does not match with BGF file.")

    # 4. Update coordinates from the snapshot
    for t in tqdm.tqdm(timesteps, ncols=120, desc="Calculating graphene z positions"):
        mybgf = get_timestep(mybgf, trj_file, t)
        avg_gra_z1 = bt.atoms_average(mybgf, 'atom.z', selection="'C_2G' in atom.ffType and atom.rNo == 1")
        avg_gra_z2 = bt.atoms_average(mybgf, 'atom.z', selection="'C_2G' in atom.ffType and atom.rNo == 2")
        dist = avg_gra_z2 - avg_gra_z1
        eff_dist = avg_gra_z2 - avg_gra_z1 - 2 * r_vdw_C
        result[t] = [avg_gra_z1, avg_gra_z2, dist, eff_dist]

    # 5. Analyze
    timesteps = sorted(result.keys())
    with open('z_position.test.dat', 'w') as f:
        output = ''
        for t in timesteps:
            output += "%d " % t
            output += "%8.3f %8.3f %8.3f %8.3f\n" % result[t]

        f.write(output)

    # 7. Done


if __name__ == "__main__":
    bgf_file = ""; trj_file = ""; ff_file = ""; z_ranges = '';

    usage = "%s -b bgf_file -t trj_file -f ff_file" % sys.argv[0]
    if len(sys.argv) < 2:
        print(usage)
        sys.exit(0)

    options, args = getopt.getopt(sys.argv[1:], 'hb:t:f:', ['help', 'bgf=', 'trj=', 'ff='])

    print(sys.argv)
    print("Requested options: " + str(options))

    for option, value in options:
        if option in ('-h', '--help'):
            print(usage)
            sys.exit(0)
        elif option in ('-b', '--bgf'):
            bgf_file = value
        elif option in ('-t', '--trj'):
            trj_file = value
        elif option in ('-f', '--ff'):
            ff_file = str(value).strip()
        elif option == NULL:
            print(usage)
            sys.exit(0)

    count_layer(bgf_file, trj_file, ff_file=ff_file)
