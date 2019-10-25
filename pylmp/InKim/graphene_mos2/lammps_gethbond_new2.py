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

def getHbond(bgf_file, trj_file, ff_file='', selection='', out_file=''):
    '''analyze something within a timestep in a series of lammps trajectory.
    '''
    # variables
    result = dict()

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
    timesteps = mytrj.load()
    N_HEADER = mytrj.nheader
    N_ATOMS = mytrj.natoms[timesteps[0]]
    N_BUFFER = N_HEADER + N_ATOMS
    if N_BGF_ATOMS != N_ATOMS:
        nu.die("Number of atoms in trajectory file does not match with BGF file.")

    # 3. Determine dump style
    dump_keywords = mytrj.dumpstyle
    yes_scale = False
    if 'xs' in dump_keywords:
        yes_scale = True

    # 4. Update coordinates from the snapshot
    dump = get_line(trj_file)
    for t in tqdm.tqdm(timesteps, ncols=120, desc="Analyzing HBonds"):
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
        def calc_hbonds():
            # variables
            A = []; D = []; hbonds = []; 
            d_crit = 3.5; a_crit = 30.0

            for atom in mybgf.a:
                if selection:
                    if "O" in atom.ffType and eval(selection):
                        A.append(atom)
                    if "O" in atom.ffType and eval(selection):
                        D.append(atom)
                else:
                    if "O" in atom.ffType:
                        A.append(atom)
                    if "O" in atom.ffType:
                        D.append(atom)
            if not len(A) or not len(D):
                nu.die("There are no atoms which can make H_bond (O atoms so far)!")

            # calculate hbonds
            for d_atom in D:
                d = np.array([d_atom.x, d_atom.y, d_atom.z])
                neigh_anos = bt.get_neighbors_aNo(A, d, r=d_crit, pbc=mytrj.pbc[t])
                for ano in neigh_anos:
                    a_atom = mybgf.getAtom(ano)
                    h = bt.is_hbonded2(mybgf, d_atom, a_atom)  # returns hbonded H coord
                    if len(h) != 0:
                        a = np.array([a_atom.x, a_atom.y, a_atom.z])
                        dist = nu.pbc_dist(a, d, mytrj.pbc[t])
                        u = h - d; v = a - d; 
                        theta = np.dot(u, v) / norm(u) / norm(v); theta = np.degrees(arccos(theta))
                        hbonds.append([d_atom.aNo, a_atom.aNo, d, a, dist, theta])

            return hbonds

        hbonds = calc_hbonds()
        result[t] = hbonds

        #break; # tester

    # 5. Analyze
    if not out_file:
        out_file = trj_file + ".hbonds.new2"
    pickle_file = out_file + ".pickle"
    with open(pickle_file, 'wb') as f:
        pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)
        print("Success to save the result to a pickle file %s" % pickle_file)

    # 7. Done


if __name__ == "__main__":
    bgf_file = ""; trj_file = ""; ff_file = ""; condition = ''; out_file = ''

    usage = "%s -b bgf_file -t trj_file -f ff_file -c selection -o out_file" % sys.argv[0]
    if len(sys.argv) < 2:
        print(usage)
        sys.exit(0)

    options, args = getopt.getopt(sys.argv[1:], 'hb:t:c:f:o:', ['help', 'bgf=', 'trj=', 'condition=', 'ff=', 'out='])

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
        elif option in ('-c', '--selection'):
            condition = str(value)
        elif option in ('-o', '--out'):
            out_file = str(value)
        elif option == NULL:
            print(usage)
            sys.exit(0)

    getHbond(bgf_file, trj_file, selection=condition, out_file=out_file)
