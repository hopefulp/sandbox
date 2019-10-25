#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys
import getopt
import optparse
import cPickle as pickle
import itertools

import numpy as np
from numpy import arccos, arcsin
from numpy.linalg import norm
import tqdm

import bgf
import bgftools as bt
import lammpstrj as lt
import nutils as nu

# Constants
O_epsilon = 0.15530000;    # SPC/E OW in kcal/mol
O_sigma = 3.16598677;      # SPC/E OW in A

def update_coord(chunk, mybgf, pbc="", scaled=False):

    coords = chunk[9:]

    for c in coords:
        c = c.split(' ')
        atom = mybgf.getAtom(int(c[0]))

        if scaled:
            atom.x = float(c[2]) * pbc[0]
            atom.y = float(c[3]) * pbc[1]
            atom.z = float(c[4]) * pbc[2]
        else:
            atom.x = float(c[2])
            atom.y = float(c[3])
            atom.z = float(c[4])

    return mybgf


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
    atom_frags = bt.getMoleculeList(mybgf)

    # 2. Read LAMMPS Trajectory
    mytrj = lt.lammpstrj(trj_file)
    mytrj.load()
    timesteps = sorted(mytrj.timesteps)
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
    for t in tqdm.tqdm(timesteps, ncols=120, desc="Analyzing HBonds"):
        chunk = [next(dump) for i in range(N_BUFFER)]

        mybgf = update_coord(chunk, mybgf, mytrj.pbc[t], scaled=yes_scale)
        mybgf = bt.periodicMoleculeSort(mybgf, mybgf.CRYSTX, fragments=atom_frags, ff_file=ff_file, silent=True)

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
                d = np.array([d_atom.x, d_atom.y, d_atom.z])    # donor coord
                neigh_anos = bt.get_neighbors_aNo(A, d, r=d_crit, pbc=mytrj.pbc[t], k=6)
                donors = [d_atom.aNo] + d_atom.CONECT

                for ano in neigh_anos:
                    a_atom = mybgf.getAtom(ano)
                    a = np.array([a_atom.x, a_atom.y, a_atom.z])    # acceptor coord
                    acceptors = [a_atom.aNo] + a_atom.CONECT

                    for ano in d_atom.CONECT:
                        h_atom = mybgf.getAtom(ano)
                        h = np.array([h_atom.x, h_atom.y, h_atom.z])
                        u = h - d; v = a - d; 
                        theta = np.dot(u, v) / norm(u) / norm(v); theta = np.degrees(arccos(theta))
                        if theta < a_crit:  # HBond exists
                            dist = nu.pbc_dist(a, d, mytrj.pbc[t])
                            dist_ah = nu.pbc_dist(d, h, mytrj.pbc[t])

                            # E_vdw
                            sigma_r = O_sigma / dist; sigma_r_6 = sigma_r**6; sigma_r_12 = sigma_r**12
                            E_vdw = 4.0 * O_epsilon * (sigma_r_12 - sigma_r_6); # E_vdw in kcal/mol

                            # E_coul
                            E_coul = 0.0
                            for i, j in itertools.product(donors, acceptors):
                                atom1 = mybgf.getAtom(i)
                                atom2 = mybgf.getAtom(j)
                                a1 = [atom1.x, atom1.y, atom1.z]
                                a2 = [atom2.x, atom2.y, atom2.z]
                                dist_ij = nu.pbc_dist(a1, a2, mytrj.pbc[t])
                                E_coul += 332.06371 * atom1.charge * atom2.charge / dist_ij # E_coul in kcal/mol

                            # E_hbond
                            E_hbond = E_coul + E_vdw  # E_hbond = E_vdw + E_coul

                            # update for v4
                            # angle between H-O-H plane and O..O vector
                            H1 = mybgf.getAtom(d_atom.CONECT[0]); h1 = [H1.x, H1.y, H1.z]    # H1
                            H2 = mybgf.getAtom(d_atom.CONECT[1]); h2 = [H2.x, H2.y, H2.z]    # H2
                            p = d - h1; q = d - h2; n = np.cross(p, q)  # normal vector of the H1-O-H2 plane
                            m = a - d;  # O..O vector
                            alpha = np.dot(n, m) / norm(n) / norm(m); alpha = np.degrees(arcsin(alpha)) # angle between H-O-H plane and O..O vector

                            #hbonds.append([d_atom.aNo, a_atom.aNo, d, a, dist, theta, [E_coul, E_vdw, E_hbond]])  # v2
                            #hbonds.append([d_atom.aNo, a_atom.aNo, d, a, dist, theta, [E_coul, E_vdw, E_hbond], dist_ah])   # v3
                            hbonds.append([d_atom.aNo, a_atom.aNo, d, a, dist, theta, [E_coul, E_vdw, E_hbond], dist_ah, alpha])   # v4

            return hbonds

        hbonds = calc_hbonds()
        result[t] = hbonds

        #break; # tester

    # 5. Analyze
    if not out_file:
        out_file = trj_file + ".hbonds.v4"
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
