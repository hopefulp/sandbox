"""
A template to write an analysis code iterating over a LAMMPS trajectory upon a BGF file.

The code is written in Python 3 and supports a multi-threaded CPU.
In python 2, this code might not work due to incomplete support of multi-threading.

Usage: modify run() with your own analysis code.
"""
# Python 3 basic packages
import sys
import _pickle as pickle
import concurrent.futures
import multiprocessing
import getopt
import logging

# additional packages with pip
import tqdm

# custom codes
import bgf
import bgftools as bt
from lammps.trj import Trj

version = "20170214"    # updated

# workers
try:
    workers = multiprocessing.cpu_count()
except NotImplementedError:
    workers = 1

def get_line(file):
    """
    This method is used to open a lengthy file as a lazy stream.
    :param file: a string for a text filename
    :return: yield a line
    """
    assert isinstance(file, str), "Wrong filename specified."
    with open(file, 'r') as myf:
        for line in myf:
            yield line


def update_coord(chunk, system, pbc, scaled=False):
    """
    This method is used to update coordinates of a BGF object from a trajectory.
    :param chunk: a list of coordinates from get_line
    :param system: a bgf.BgfFile class to update coordinates
    :param pbc: a list [x, y, z] of pbc parameters
    :param scaled: yes if coordinates are from scaled LAMMPS trajectory (i.e. xs, ys, and zs)
    :return: a bgf.BgfFile class (with updated coordinates)
    """
    assert isinstance(system, bgf.BgfFile), "Wrong object: should be an instance of bgf.BgfFile class."
    assert len(pbc) == 3, "PBC list should have length 3."

    # write pbc
    system.CRYSTX[0] = pbc[0]
    system.CRYSTX[1] = pbc[1]
    system.CRYSTX[2] = pbc[2]

    coords = chunk[9:]

    for c in coords:
        c = c.split(' ')
        atom = system.getAtom(int(c[0]))

        if scaled:
            atom.x = float(c[2]) * pbc[0]
            atom.y = float(c[3]) * pbc[1]
            atom.z = float(c[4]) * pbc[2]
        else:
            atom.x = float(c[2])
            atom.y = float(c[3])
            atom.z = float(c[4])

    return system


def run(system, t, coords, sel=""):
    """
    This function performs the analysis.
    Input:
        - bgf.BgfFile mybgf
        - int t: a timestep written in LAMMPS trjectory
        - list chunk: list of atom coordinates in LAMMPS trajectory
        - str sel: a string to evaluate for selecting atoms
    Output:
        - int t: timestep
        - list result: found hbonds
    """
    assert isinstance(system, bgf.BgfFile), "Wrong object: should be an instance of bgf.BgfFile class."
    assert isinstance(t, int), "Timestep should be int."
    assert isinstance(coords, list), "coords should be a list."

    # update the coordinates from trj to the BGF object
    system = update_coord(coords, system, mytrj.pbc[t], scaled=yes_scale)

    # wrap the coordinates into the pbc cell
    system = bt.periodicMoleculeSort(system, system.CRYSTX, fragments=atom_frags, ff_file=ff_file, silent=True)

    selected_atoms = [atom for atom in system.a if eval(sel)]

    # do analysis here
    out = []  # list to save results
    # ...
    # something something something
    # ...
    #
    # blah blah blah
    # ...

    # save each trajectory as BGF file if verbose mode activated.
    if verbose:
        system.saveBGF(str(t) + ".bgf")

    return t, out


if __name__ == "__main__":

    bgf_file = ""
    trj_file = ""
    ff_file = ""
    selection = ''
    out_file = ''
    verbose = False
    n_sample = 0
    n_process = 0

    usage = "%s -b bgf_file -t trj_file -f ff_file -s selection -o out_file -n n_timesteps" % sys.argv[0]

    # logger
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s: %(message)s", datefmt="%y-%m-%d %H:%M")
    logger = logging.getLogger("Analyzer")

    if len(sys.argv) < 2:
        logger.info(usage)
        sys.exit(0)

    options, args = getopt.getopt(sys.argv[1:], 'hb:t:s:f:o:vn:p:',
                                  ['help', 'bgf=', 'trj=', 'selection=', 'ff=', 'out=', 'verbose=', 'n=', 'process='])

    logger.info(sys.argv)
    logger.info("Requested options: " + str(options))
    logger.info("Initializing multiprocessing: Found %d processors" % workers)

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
            logger.setLevel(logging.DEBUG)
        elif option in ('-n', '--n'):
            n_sample = int(value)
        elif option in ('-p', '--process'):
            n_process = int(value)
        else:
            print(usage)
            sys.exit(0)

    # mandatory parameters
    assert bgf_file, "No BGF file specified."
    assert trj_file, "No LAMMPS trajectory file specified."

    # options
    if not out_file:
        out_file = "myanalyze"

    # 1. Load BGF
    mybgf = bgf.BgfFile(bgf_file)
    N_BGF_ATOMS = len(mybgf.a)
    atom_frags = bt.getMoleculeList(mybgf)

    # 2. Read LAMMPS Trajectory
    mytrj = Trj(trj_file)
    mytrj.load()
    timesteps = sorted(mytrj.timesteps)
    N_HEADER = mytrj.nheader
    N_ATOMS = mytrj.natoms[timesteps[0]]
    N_BUFFER = N_HEADER + N_ATOMS
    if N_BGF_ATOMS != N_ATOMS:
        logger.error("Number of atoms in trajectory file does not match with BGF file.")
        sys.exit(0)

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
        logger.info("Only last %d steps will be calculated." % n_sample)
        timesteps = timesteps[-n_sample:]

    # 5. Analyze
    result = []
    if workers == 1:    # single processor usage
        logger.info("\nStarting job.. Will use %d processors for calculation." % workers)
        for t in tqdm.tqdm(timesteps, ncols=120):
            result = run(mybgf, t, chunks[t])
    else:   # multi processors usage
        if workers > 3:
            actual_workers = workers - 1
        else:
            actual_workers = workers

        if n_process:  # number of processor specified
            actual_workers = n_process

        logger.info("\nStarting job.. Will use %d processors in parallel for calculation." % actual_workers)
        with concurrent.futures.ProcessPoolExecutor(max_workers=actual_workers) as exe:
            fs = {exe.submit(run, mybgf, t, chunks[t], selection) for t in timesteps}
            kwargs = {'total': len(fs), 'unit': 'run', 'leave': True, 'ncols': 120, 'desc': "Analyzing snapshots"}

            for f in tqdm.tqdm(concurrent.futures.as_completed(fs), **kwargs):
                pass

            done, _ = concurrent.futures.wait(fs)
            result = [f.result() for f in done]

    result = {r[0]: r[1] for r in result}

    # 6. Save results
    pickle_file = out_file + ".analysis.pickle"

    with open(pickle_file, 'wb') as f:
        logger.info("Writing results to the pickle file %s" % pickle_file)
        pickle.dump(result, f)
        logger.info("Success to save the result to a pickle file %s" % pickle_file)

    logger.info("Done.")
    sys.exit(0)