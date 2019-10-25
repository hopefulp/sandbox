#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys
import getopt

import tqdm

import bgf
import bgftools as bt
import lammpstrj as lt
import nutils as nu
"""
REMARK: Benchmark shows that reading file in this method gradually increases the cost to fetch coordinates from the trajectory file as t proceeds.
"""

def get_timestep(mybgf, trj_file, request_t, ff_file='', out_file='', wrap=False):
    '''get a set of coordinates of atoms for a requested timestep request_t.
    '''

    # 2. Read LAMMPS Trajectory
    mytrj = lt.lammpstrj(trj_file)
    timesteps = mytrj.load()
    N_BGF_ATOMS = len(mybgf.a)
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

    def get_line(file):
        with open(file, 'r') as f:
            for line in f:
                yield line

    # 4. Update coordinates from the snapshot
    dump = get_line(trj_file)
    for t in timesteps:
        if t == request_t:
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

            if wrap:
                mybgf = bt.periodicMoleculeSort(mybgf, mybgf.CRYSTX, ff_file=ff_file, silent=True)

            return mybgf

        else:
            # skip reading
            for i in range(N_BUFFER):
                next(dump)

    nu.warn("Could not load the trajectory.")
    return False
