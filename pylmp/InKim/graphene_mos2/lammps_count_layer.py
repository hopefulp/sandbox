#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys
import getopt
import optparse
if sys.version_info < (3, 0):
    import cPickle as pickle
else:
    import _pickle as pickle

import numpy as np
from numpy import arccos
from numpy.linalg import norm
import tqdm

import bgf
import bgftools as bt
import lammpstrj as lt
import nutils as nu

def count_layer(bgf_file, trj_file, z_range, ff_file='', out_file='', selection=''):
    '''counts number of O atoms in each layer.
    '''
    # variables
    result = dict()
    z_bin = [float(i) for i in z_range.split()]

    # inner functions
    def get_line(file):
        with open(file, 'r') as f:
            for line in f:
                yield line

    # 1. Load BGF
    mybgf = bgf.BgfFile(bgf_file)
    N_BGF_ATOMS = len(mybgf.a)

    # 2. Read LAMMPS Trajectory
    #timesteps, N_HEADER, N_ATOMS = lt.getTrjInfo(trj_file)
    mytrj = lt.lammpstrj(trj_file)
    timesteps = mytrj.timesteps
    N_HEADER = mytrj.nheader
    N_ATOMS = mytrj.natoms[timesteps[0]]
    N_BUFFER = N_HEADER + N_ATOMS
    if N_BGF_ATOMS != N_ATOMS:
        nu.die("Number of atoms in trajectory file does not match with BGF file.")

    # 3. Determine dump style
    dump_keywords = mytrj._dump_style
    yes_scale = False
    if 'xs' in dump_keywords:
        yes_scale = True

    # 4. Update coordinates from the snapshot
    dump = get_line(trj_file)
    for t in tqdm.tqdm(timesteps, ncols=120, desc="Counting Atoms"):
        chunk = [next(dump) for i in range(N_BUFFER)]

        t = int(chunk[1])
        mybgf.CRYSTX = mytrj.pbc[t] + [90.0, 90.0, 90.0]

        coords = chunk[9:]
        for c in coords:
            c = c.split(' ')
            atom = mybgf.getAtom(int(c[0]))

            if yes_scale:
                atom.x = float(c[2]) * mytrj.pbc[t][0]
                atom.y = float(c[3]) * mytrj.pbc[t][1]
                atom.z = float(c[4]) * mytrj.pbc[t][2]
            else:
                atom.x = float(c[2])
                atom.y = float(c[3])
                atom.z = float(c[4])

        #mybgf = bt.periodicMoleculeSort(mybgf, mybgf.CRYSTX, ff_file=ff_file, silent=True)

        ### collect data
        def layer():
            # variables
            A = []; layer = dict()
            if selection:
                A = [atom for atom in mybgf.a if ("O" in atom.ffType and atom.chain == "I" and eval(selection))]
            else:
                A = [atom for atom in mybgf.a if "O" in atom.ffType and atom.chain == "I"]

            def find_layer(z):
                if z_bin:
                    return int(np.digitize(z, z_bin))
                else:
                    return 0

            # find layer
            for atom in A:
                layer[atom.aNo] = find_layer(atom.z)

            return layer

        layer = layer()
        result[t] = layer

    # 5. Analyze
    pickle_file = "layer." + z_range.replace(' ', '_').replace('.', 'p') + ".pickle"
    pickle_file = pickle_file.replace('..', '.')
    with open(pickle_file, 'wb') as f:
        pickle.dump(result, f)
        print("Success to save the result to a pickle file %s" % pickle_file)

    # 7. Done


if __name__ == "__main__":
    bgf_file = ""; trj_file = ""; ff_file = ""; z_ranges = ''; sel = ""

    usage = "%s -b bgf_file -t trj_file -f ff_file -r 'z_bin_ranges' -s 'selection'" % sys.argv[0]
    if len(sys.argv) < 2:
        print(usage)
        sys.exit(0)

    options, args = getopt.getopt(sys.argv[1:], 'hb:t:r:f:s:', ['help', 'bgf=', 'trj=', 'ranges=', 'ff=', 'sel='])

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
        elif option in ('-r', '--ranges'):
            z_ranges = str(value)
        elif option in ('-s', '--sel'):
            sel = value
        elif option == NULL:
            print(usage)
            sys.exit(0)

    count_layer(bgf_file, trj_file, z_ranges, ff_file=ff_file, selection=sel)
