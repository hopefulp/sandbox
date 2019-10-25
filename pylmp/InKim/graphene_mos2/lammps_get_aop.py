#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys
import getopt
import cPickle as pickle
import itertools

import numpy as np
import numpy.linalg as la
import scipy.spatial

import periodic_kdtree as pkdtree
import tqdm

import bgf
import bgftools as bt
import lammpstrj as lt
import nutils as nu

version = '20160527'

def analyze(bgf_file, trj_file, ff_file='', name='', selection='', d_crit=3.5):
    '''analyze something within a timestep in a series of lammps trajectory.
    '''
    # variables
    result_angles = dict()
    result_diheds = dict()
    print("Using neighboring water distance criteria %4.1f A (should be consistent with the coordination number)" % d_crit)

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
    N_ATOMS = mytrj.natoms[0]
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
    for t in tqdm.tqdm(timesteps, ncols=120, desc="Analyzing"):
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

        # Calculate Angles and Dihedrals for Water Structure
        O = [];
        for atom in mybgf.a:
            if selection:
                if "O" in atom.ffType and eval(selection):
                    O.append(atom)
            else:
                if "O" in atom.ffType:
                    O.append(atom)

        if not len(O):
            nu.warn("There are no O atoms which satisfies %s!" % selection)

        coords = []; anos = [];
        for atom in O:
            coords.append([atom.x, atom.y, atom.z])
            anos.append(atom.aNo)

        coords = np.array(coords)
        #tree = pkdtree.PeriodicKDTree(mytrj.pbc[t], coords)    # KDTree with pbc
        tree = scipy.spatial.KDTree(coords, leafsize=len(coords)+1)   # REMARK: normal KDTree -- do not consider over pbc

        angles = []; diheds = [];
        for atom in O:
            # Find neighbors
            coord = np.array([atom.x, atom.y, atom.z])
            neighbor_O = []; # list of O atoms near centerO
            d, ndx = tree.query(coord, k=5)
            index = np.where((d >= 1e-4) & (d <= d_crit))  # discard self
            for i in ndx[index]:
                neighbor_O.append(mybgf.getAtom(anos[i]))

            # Angles
            for i, j in itertools.combinations(neighbor_O, 2):
                x = [i.x, i.y, i.z]; y = [j.x, j.y, j.z]
                angle = nu.angle(x, coord, y, radians=False)
                angles.append([angle, i.z, atom.z, j.z])    # angles: result_angle data structure


            # Dihedrals: farthest H-O...O-H for all neighboring O-O pairs
            for atom2 in neighbor_O:
                # are they have an hydrogen bond?
                if bt.is_hbonded(mybgf, atom, atom2):
                    hO_anos = atom.CONECT   # H connected to center O
                    hO2_anos = atom2.CONECT # H connected to neighboring O
                    max_dist_pair = []; max_dist = 0.0;
                    for k, l in itertools.product(hO_anos, hO2_anos):
                        h1 = mybgf.getAtom(k)
                        h2 = mybgf.getAtom(l)
                        dist = nu.dist([h1.x, h1.y, h1.z], [h2.x, h2.y, h2.z])
                        if dist > max_dist:
                            max_dist = dist
                            max_dist_pair = [h1, h2]

                phi = bgf.dihedral(max_dist_pair[0], atom, atom2, max_dist_pair[1], radians=False)
                diheds.append([phi, atom.z, atom2.z])   # diheds: result_dihed data structure

        result_angles[t] = angles
        result_diheds[t] = diheds

    # Write-ups
    if name:
        out_file = name + '.aop'
    else:
        out_file = 'aop'

    angle_out_file = out_file + ".angle.pickle"
    dihed_out_file = out_file + ".dihed.pickle"

    with open(angle_out_file, 'wb') as f:
        pickle.dump(result_angles, f, protocol=pickle.HIGHEST_PROTOCOL)
        print("Success to save the angle result to a pickle file %s" % angle_out_file)

    with open(dihed_out_file, 'wb') as f:
        pickle.dump(result_diheds, f, protocol=pickle.HIGHEST_PROTOCOL)
        print("Success to save the dihed result to a pickle file %s" % dihed_out_file)

    print("Done.")
    
    # Done


if __name__ == "__main__":
    bgf_file = ""; trj_file = ""; ff_file = ""; task_name = ""; distance = 0.0; sel = "";

    usage = "%s -b bgf_file -t trj_file -f ff_file -s task_name" % sys.argv[0]
    if len(sys.argv) < 2:
        print(usage)
        sys.exit(0)

    options, args = getopt.getopt(sys.argv[1:], 'hb:t:f:s:d:c:', ['help', 'bgf=', 'trj=', 'ff=', 'name=', 'distance=', 'selection='])

    print "Requested options: " + str(options)

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
        elif option in ('-s', '--name'):
            prefix = value
        elif option in ('-d', '--distance'):
            distance = float(value)
        elif option in ('-c', '--selection'):
            sel = value
        elif option == NULL:
            print(usage)
            sys.exit(0)

    if not distance:
        distance = 3.5

    analyze(bgf_file, trj_file, name=task_name, ff_file=ff_file, d_crit=distance, selection=sel)
