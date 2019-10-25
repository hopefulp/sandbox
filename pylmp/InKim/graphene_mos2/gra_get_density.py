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
import nutils as nu

def count_layer(bgf_file, trj_file, ff_file='', out_file='', prefix=''):
    '''measures confined density of graphene interlayer
    '''
    # variables
    result = dict()
    r_vdw_C = 3.38383824/2

    # inner functions
    def get_line(file):
        with open(file, 'r') as f:
            for line in f:
                yield line

    # 1. Load BGF
    mybgf = bgf.BgfFile(bgf_file)
    N_BGF_ATOMS = len(mybgf.a)

    # 2. Read LAMMPS Trajectory
    mytrj = lt.lammpstrj(trj_file)
    timesteps = mytrj.load()
    N_HEADER = mytrj.nheader
    N_ATOMS = mytrj.natoms[0]
    N_BUFFER = N_HEADER + N_ATOMS
    if N_BGF_ATOMS != N_ATOMS:
        nu.die("Number of atoms in trajectory file does not match with BGF file.")

    # 3. Determine dump style
    dump_keywords = mytrj.dumpstyle
    yes_scale = False
    if 'xs' in dump_keywords:
        yes_scale = True

    # Mass
    wat_ano = []
    for atom in mybgf.a:
        if atom.chain == "I":
            wat_ano.append(atom.aNo)

    mass = bt.getMass(mybgf, wat_ano, ff_file=ff_file)

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
                atom.x = float(c[2]) * pbc[0]
                atom.y = float(c[3]) * pbc[1]
                atom.z = float(c[4]) * pbc[2]
            else:
                atom.x = float(c[2])
                atom.y = float(c[3])
                atom.z = float(c[4])

        mybgf = bt.periodicMoleculeSort(mybgf, mybgf.CRYSTX, ff_file=ff_file, silent=True)

        ### collect data
        def density():
            # variables
            density = dict()

            # height 
            avg_gra_z1 = bt.atoms_average(mybgf, 'atom.z', selection="'C_2G' in atom.ffType and atom.rNo == 1")
            avg_gra_z2 = bt.atoms_average(mybgf, 'atom.z', selection="'C_2G' in atom.ffType and atom.rNo == 2")
            dist = avg_gra_z2 - avg_gra_z1
            eff_dist = avg_gra_z2 - avg_gra_z1 - 2 * r_vdw_C

            # volume
            x = mytrj.pbc[t][0]
            y = mytrj.pbc[t][1]
            volume = x * y * dist
            eff_volume = x * y * eff_dist

            # density
            density = mass / 6.022 / volume * 10
            eff_density = mass / 6.022 / eff_volume * 10

            return [x, y, dist, mass, volume, eff_volume, eff_dist, density, eff_density]

        density = density()
        result[t] = density

    # 5. Analyze
    timesteps = sorted(result.keys())
    with open(prefix + '.trj.density.dat', 'w') as f:
        output = "#x\ty\tdist\tmass\tvolume\teff_volume\teff_dist\tdensity\teff_density\n"
        for t in timesteps:
            output += "%d " % t
            output += " ".join("%8.3f" % float(i) for i in result[t])

        f.write(output)


    # 7. Done


if __name__ == "__main__":
    bgf_file = ""; trj_file = ""; ff_file = ""; prefix = '';

    usage = "%s -b bgf_file -t trj_file -f ff_file -p prefix" % sys.argv[0]
    if len(sys.argv) < 2:
        print(usage)
        sys.exit(0)

    options, args = getopt.getopt(sys.argv[1:], 'hb:t:p:f:', ['help', 'bgf=', 'trj=', 'prefix=', 'ff='])

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
        elif option in ('-p', '--prefix'):
            prefix = value
        elif option == NULL:
            print(usage)
            sys.exit(0)

    count_layer(bgf_file, trj_file, ff_file=ff_file, prefix='')
