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
    print("Found %d processors" % workers)

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
        A = [atom for atom in mybgf.a if "O" in atom.ffType and eval(selection)]
        D = [atom for atom in mybgf.a if "O" in atom.ffType and eval(selection)]
    else:
        A = [atom for atom in mybgf.a if "O" in atom.ffType]
        D = [atom for atom in mybgf.a if "O" in atom.ffType]

    if not len(A) or not len(D):
        nu.warn("There are no atoms which can make hbonds (O atoms so far)!")
        return

    # calculate hbonds
    hbonds = []; 

    for d_atom in D:
        d = np.array([d_atom.x, d_atom.y, d_atom.z])    # donor coord
        neigh_anos = bt.get_neighbors_aNo(A, d, r=d_crit, pbc=pbc, k=5)
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
                    # calculate pair energy as Hbond energy
                    dist = nu.pbc_dist(a, d, pbc)
                    dist_ah = nu.pbc_dist(d, h, pbc)

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
                        dist_ij = nu.pbc_dist(a1, a2, pbc)
                        E_coul += 332.06371 * atom1.charge * atom2.charge / dist_ij # E_coul in kcal/mol

                    # E_hbond
                    E_hbond = E_coul + E_vdw  # E_hbond = E_vdw + E_coul

                    '''
                    # calculate dreiding style hbond energy:  http://lammps.sandia.gov/doc/pair_hbond_dreiding.html
                    # REMARK: tested on bulk water, ice, and confined water between MoS2 but achieved suspicious values:
                    # 6.0 $\AA$ maximum: -7.925 kcal/mol
                    # 8.0 $\AA$ maximum: -7.875 kcal/mol
                    # 11.0 $\AA$ maximum: -8.025 kcal/mol
                    # ice maximum: -8.325 kcal/mol
                    # water maximum: -7.775 kcal/mol
                    sigma = 2.75 # in A
                    epsilon = 9.0 # in kcal/mol
                    u = d - h; v = a - h; cos_theta = np.dot(u, v) / norm(u) / norm(v); #cos_theta = np.cos(theta)
                    sigma_r = sigma / dist; sigma_r_2 = sigma_r**2; sigma_r_12 = sigma_r_2**6; sigma_r_10 = sigma_r_2**5
                    E_hbond = epsilon * ( 5 * sigma_r_12 - 6 * sigma_r_10 ) * cos_theta**4
                    '''

                    hbonds.append([d_atom.aNo, a_atom.aNo, dist, theta, E_hbond])   # v5

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

    if not sel:
        sel = calculate_sel_kwac(mybgf)    # for Kwac's confined water molecules in MoS2 nanoslit
    elif sel == "all":
        sel = ""    # for water

    if sel:
        sel_O_atoms = [atom.aNo for atom in mybgf.a if eval(sel) and 'O' in atom.ffType]
    else:
        sel_O_atoms = [atom.aNo for atom in mybgf.a if 'O' in atom.ffType]

    result = hbonds(mybgf, selection=sel)

    if verbose:
        mybgf.saveBGF(str(t) + ".bgf")

    return t, sel_O_atoms, result

    
if __name__ == "__main__":

    bgf_file = ""; trj_file = ""; ff_file = ""; selection = ''; out_file = ''; verbose = False; nsample = 0;

    usage = "%s -b bgf_file -t trj_file -f ff_file -s selection -o out_file" % sys.argv[0]
    if len(sys.argv) < 2:
        print(usage)
        sys.exit(0)

    options, args = getopt.getopt(sys.argv[1:], 'hb:t:s:f:o:vn:', ['help', 'bgf=', 'trj=', 'selection=', 'ff=', 'out=', 'verbose=', 'n='])

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
            nsample = int(value)
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
    for t in tqdm.tqdm(timesteps, ncols=120, desc="Updating coordinates"):
        chunk = [next(dump) for i in range(N_BUFFER)]
        chunks[t] = chunk

    # 4. Determine the number of trajectories specified
    if nsample:
        print("Only last %d steps will be calculated." % nsample)
        timesteps = timesteps[-nsample:]

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
        # r = (t, result, mybgf, sel)
        result2[r[0]] = r[1:]
    result = result2

    # 6. Write file
    pickle_file = out_file + ".hbond.v5.pickle"

    with open(pickle_file, 'wb') as f:
        print("Writing results to the pickle file %s" % pickle_file)
        pickle.dump(result, f)
        print("Success to save the result to a pickle file %s" % pickle_file)

    print("Done.")
