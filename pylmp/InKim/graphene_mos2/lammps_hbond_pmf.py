#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys
import getopt
import optparse
import cPickle as pickle
import itertools

import numpy as np
from numpy import arccos
from numpy.linalg import norm
import tqdm

import bgf
import bgftools as bt
import lammpstrj as lt
import nutils as nu

def getHbond(bgf_file, trj_file, ff_file='', selection='', out_file='', nsamples=0):
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
    if nsamples: 
        requested_timesteps = timesteps[-nsamples:]
    else:
        requested_timesteps = timesteps
    N_HEADER = mytrj.nheader
    N_ATOMS = mytrj.natoms[timesteps[0]]
    N_BUFFER = N_HEADER + N_ATOMS
    if N_BGF_ATOMS != N_ATOMS:
        nu.die("Number of atoms in trajectory file does not match with BGF file.")

    A = []; D = [];
    if not selection:
        for atom in mybgf.a:
            if "O" in atom.ffType:
                A.append(atom)
            if "O" in atom.ffType:
                D.append(atom)

    # 3. Determine dump style
    dump_keywords = mytrj.dumpstyle

    # 4. Update coordinates from the snapshot
    dump = get_line(trj_file)
    for t in tqdm.tqdm(timesteps, ncols=120, desc="Analyzing HBonds"):
        chunk = [next(dump) for i in range(N_BUFFER)]

        t = int(chunk[1])
        if not t in requested_timesteps:
            continue;

        mybgf.CRYSTX = mytrj.pbc[t] + [90.0, 90.0, 90.0]

        coords = chunk[9:]
        for c in coords:
            c = c.split(' ')
            atom = mybgf.getAtom(int(c[0]))

            atom.x = float(c[2])
            atom.y = float(c[3])
            atom.z = float(c[4])

        mybgf = bt.periodicMoleculeSort(mybgf, mybgf.CRYSTX, ff_file=ff_file, silent=True)

        ### collect data
        def calc_hbonds():
            hbonds = []; 
            d_crit = 3.6;

            if selection:
                for atom in mybgf.a:
                    if "O" in atom.ffType and eval(selection):
                        A.append(atom)
                    if "O" in atom.ffType and eval(selection):
                        D.append(atom)

            # calculate hbonds
            for d_atom in D:
                d = np.array([d_atom.x, d_atom.y, d_atom.z])    # donor coord
                neigh_anos = bt.get_neighbors_aNo(A, d, r=d_crit, pbc=mytrj.pbc[t], k=6)

                for ano in neigh_anos:
                    a_atom = mybgf.getAtom(ano)
                    a = np.array([a_atom.x, a_atom.y, a_atom.z])    # acceptor coord
                    dist = nu.pbc_dist(a, d, mytrj.pbc[t])
                    if dist > 2.0:
                        for ano in d_atom.CONECT:
                            h_atom = mybgf.getAtom(ano)
                            h = np.array([h_atom.x, h_atom.y, h_atom.z])
                            u = h - d; v = a - d; 
                            theta = np.dot(u, v) / norm(u) / norm(v); #theta = np.degrees(arccos(theta))
                            if -0.17364817766693034885171662676931 < theta < 1.0:
                                hbonds.append([dist, np.degrees(theta)])

            return hbonds

        hbonds = calc_hbonds()
        result[t] = hbonds

        #break; # tester

    # 5. Analyze
    if not out_file:
        out_file = trj_file + ".hbonds_pmf"
    pickle_file = out_file + ".pickle"
    with open(pickle_file, 'wb') as f:
        pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)
        print("Success to save the result to a pickle file %s" % pickle_file)

    # 6. Done


if __name__ == "__main__":
    bgf_file = ""; trj_file = ""; ff_file = ""; condition = ''; out_file = ''; n_samples = 0;

    usage = "%s -b bgf_file -t trj_file -f ff_file -c selection -o out_file" % sys.argv[0]
    if len(sys.argv) < 2:
        print(usage)
        sys.exit(0)

    options, args = getopt.getopt(sys.argv[1:], 'hb:t:c:f:o:n:', ['help', 'bgf=', 'trj=', 'condition=', 'ff=', 'out=', 'nsamples='])

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
        elif option in ('-n', '--nsamples'):
            n_samples = int(value)
        elif option == NULL:
            print(usage)
            sys.exit(0)

    getHbond(bgf_file, trj_file, selection=condition, out_file=out_file, nsamples = n_samples)
