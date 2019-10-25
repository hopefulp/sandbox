import sys
if sys.version_info < (3, 0):
    print("ERROR: Must be using Python 3.0 or higher to use parallel mode.")
    import cPickle as pickle
    workers = 1
else:
    import _pickle as pickle
    import concurrent.futures
    import multiprocessing

    # workers
    try:
        workers = multiprocessing.cpu_count()
    except NotImplementedError:
        workers = 1
    print("Initializing multiprocessing: Found %d processors" % workers)

import getopt
import optparse
import itertools

import numpy as np
from numpy import arccos, arcsin
from numpy.linalg import norm
import tqdm

import bgf
import bgftools as bt
import lammpstrj as lt
import nutils as nu
#import analysis

# Constants
O_epsilon = 0.15530000;    # SPC/E OW in kcal/mol
O_sigma = 3.16598677;      # SPC/E OW in A
mybgf = bgf.BgfFile() 


def update_coord(chunk, mybgf, pbc="", scaled=False):

    # write pbc
    mybgf.CRYSTX[0] = pbc[0]
    mybgf.CRYSTX[1] = pbc[1]
    mybgf.CRYSTX[2] = pbc[2]

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


def calculate_sel_kwac(mybgf):
    """
    Create a selection string by calculating some properties from mybgf to calculate hbonds of Kwac's structure
    Input:
        - bgf.BgfFile mybgf
    Output:
        - string sel
    """

    type = "gra"
    for atom in mybgf.a:
        if "Mo" in atom.ffType:
            type = "mos"

    # find inwater boundary for selection
    swr = 3.270615945 / 2
    gwr = 3.057430885 / 2
    margin = 5.0
    gwa_y = bt.atoms_average(mybgf, 'atom.y', selection="'GWA' in atom.rName")
    gwb_y = bt.atoms_average(mybgf, 'atom.y', selection="'GWB' in atom.rName")

    if type == "gra":
        avg_z1 = bt.atoms_average(mybgf, 'atom.z', selection="'C_R' in atom.ffType and 'GRA' in atom.rName")    # gra bottom
        avg_z2 = bt.atoms_average(mybgf, 'atom.z', selection="'C_R' in atom.ffType and 'GRB' in atom.rName")    # gra top
        actual_distance = avg_z2 - avg_z1
    elif type == "mos":
        avg_s3a_z2 = bt.atoms_average(mybgf, 'atom.z', selection="'S_3a' in atom.ffType and atom.rNo == 2")    # mos2 top
        avg_s3b_z1 = bt.atoms_average(mybgf, 'atom.z', selection="'S_3b' in atom.ffType and atom.rNo == 1")    # mos2 bottom
        actual_distance = avg_s3a_z2 - avg_s3b_z1

    inwater_x = mybgf.CRYSTX[0]
    inwater_y = gwb_y - gwa_y + 2 * gwr - 2 * margin
    if type == "gra":
        inwater_z = actual_distance - 2 * gwr
    elif type == "mos":
        inwater_z = actual_distance - 2 * swr

    sel = "atom.y > {gwa_y} + {margin} and atom.y < {gwb_y} - {margin}".format(**vars())

    return sel


def hbonds(mybgf, selection = ""):
    '''
    This function calculates Hydrogen bonds(hbonds) between O(donor) and O(acceptor), especially for water molecules.
    Input: 
        - bgf.BgfFile mybgf
        - string selection: a region to search donor and acceptor atoms in mybgf
    Output:
        - int n_sel_atoms: number of hbond-able atoms in the selection region
        - list hbonds: self-descriptive
    '''

    # variables
    pbc = mybgf.CRYSTX[:3]
    d_crit = 3.5; a_crit = 30.0 # Chandler's criteria

    if selection:
        D = [atom for atom in mybgf.a if ("O" in atom.ffType or "N" in atom.ffType or "F" in atom.ffType or "S" in atom.ffType) and eval(selection)]
        A = [atom for atom in mybgf.a if ("O" in atom.ffType or "N" in atom.ffType or "F" in atom.ffType or "S" in atom.ffType) and eval(selection)]
    else:
        D = [atom for atom in mybgf.a if ("O" in atom.ffType or "N" in atom.ffType or "F" in atom.ffType or "S" in atom.ffType)]
        A = [atom for atom in mybgf.a if ("O" in atom.ffType or "N" in atom.ffType or "F" in atom.ffType or "S" in atom.ffType)]

    if not len(A) or not len(D):
        nu.warn("There are no atoms which can make hbonds!")
        return

    # calculate hbonds
    hbonds = []; 

    for d_atom in tqdm.tqdm(D, ncols=120, desc='Iterating'):
        d = np.array([d_atom.x, d_atom.y, d_atom.z])    # donor coord
        donors_H = [ano for ano in d_atom.CONECT if "H" in mybgf.getAtom(ano).ffType]
        neigh_anos, distances = bt.get_neighbors_aNo(A, d, r=d_crit, pbc=pbc, k=5, return_distance=True) # atoms within d_crit of donor will be acceptor candidates

        for index, ano in enumerate(neigh_anos):
            a_atom = mybgf.getAtom(ano)
            a = np.array([a_atom.x, a_atom.y, a_atom.z])    # acceptor coord

            for ano in donors_H:
                h_atom = mybgf.getAtom(ano)
                h = np.array([h_atom.x, h_atom.y, h_atom.z])
                u = h - d; v = a - d; 
                theta = np.dot(u, v) / norm(u) / norm(v); theta = np.degrees(arccos(theta))
                if theta < a_crit:  # HBond exists
                    #dist = nu.pbc_dist(a, d, pbc)
                    #dist = distances[0][index]
                    #dist_ah = nu.pbc_dist(d, h, pbc)

                    hbonds.append([d_atom.aNo, a_atom.aNo, distances[0][index+1], theta])   # v5.2

    return hbonds


def run(mybgf, t, chunk, sel=""):
    """
    This function performs the calculation.
    Output:
        - int t: timestep
        - list result: found hbonds
        - bgf.BgfFile mybgf: BgfFile with update coordinates
        - string sel: a selection
    """

    mybgf = update_coord(chunk, mybgf, mytrj.pbc[t], scaled=yes_scale)
    mybgf = bt.periodicMoleculeSort(mybgf, mybgf.CRYSTX, fragments=atom_frags, ff_file=ff_file, silent=True)

    result = hbonds(mybgf, selection=sel)

    if verbose:
        mybgf.saveBGF(str(t) + ".bgf")

    return t, result

    
if __name__ == "__main__":

    bgf_file = ""; trj_file = ""; ff_file = ""; selection = ''; out_file = ''; verbose = False; n_sample = 0; n_process = 0;

    usage = "%s -b bgf_file -t trj_file -f ff_file -s selection -o out_file -n n_timesteps" % sys.argv[0]
    if len(sys.argv) < 2:
        print(usage)
        sys.exit(0)

    options, args = getopt.getopt(sys.argv[1:], 'hb:t:s:f:o:vn:p:', ['help', 'bgf=', 'trj=', 'selection=', 'ff=', 'out=', 'verbose=', 'n=', 'process='])

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
        elif option in ('-s', '--selection'):
            selection = str(value)
        elif option in ('-o', '--out'):
            out_file = str(value)
        elif option in ('-v', '--verbose'):
            verbose = True
        elif option in ('-n', '--n'):
            n_sample = int(value)
        elif option in ('-p', '--process'):
            n_process = int(value)
        elif option == NULL:
            print(usage)
            sys.exit(0)


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
    mytrj = lt.lammpstrj(trj_file); mytrj.load()
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
    chunks = dict()
    for t in tqdm.tqdm(timesteps, ncols=120, desc="Dumping coordinates"):
        chunk = [next(dump) for i in range(N_BUFFER)]
        chunks[t] = chunk

    # 4. Determine the number of trajectories specified
    if n_sample:
        print("Only last %d steps will be calculated." % n_sample)
        timesteps = timesteps[-n_sample:]

    # 5. Calculate hbonds
    result = []
    if workers == 1:
        # single processor usage
        print("\nStarting job.. Will use %d processors for calculation." % workers)
        for t in tqdm.tqdm(timesteps, ncols=120):
            result = run(mybgf, t, chunks[t], sel=selection)
    else:
        # multi processors usage
        if workers > 3: 
            actual_workers = workers - 1
        else:
            actual_workers = workers

        if n_process:   # number of processor specified
            actual_workers = n_process

        print("\nStarting job.. Will use %d processors in parallel for calculation." % actual_workers)
        with concurrent.futures.ProcessPoolExecutor(max_workers=actual_workers) as exe:
            #fs = {exe.submit(run, mybgf, k, chunks[k], selection, verbose) for k in chunks.keys()}
            fs = {exe.submit(run, mybgf, t, chunks[t], selection) for t in timesteps}
            kwargs = {'total': len(fs), 'unit': 'run', 'leave': True, 'ncols': 120, 'desc': "Calculating Hbonds"}

            # progress bar: from https://renzo.lucioni.xyz/parallel-progress/
            for f in tqdm.tqdm(concurrent.futures.as_completed(fs), **kwargs):
                pass

            done, _ = concurrent.futures.wait(fs)
            result = [f.result() for f in done]

    result2 = dict()
    for r in result:
        # r = (t, result)
        result2[r[0]] = r[1]
    result = result2

    # 6. Write file
    pickle_file = out_file + ".hbond.v5.2.pickle"

    with open(pickle_file, 'wb') as f:
        print("Writing results to the pickle file %s" % pickle_file)
        pickle.dump(result, f)
        print("Success to save the result to a pickle file %s" % pickle_file)

    print("Done.")
