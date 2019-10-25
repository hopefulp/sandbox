# python modules
import sys
import os
import pickle
import time
import itertools
import tqdm
import nutils as nu
from nutils import timer
import memory_footprint

version = '20160712'

class lammpstrj():

    def _get_line(self, trj_file):
        with open(trj_file, 'r') as file:
            for i in file:
                yield i

    #---------------------------

    def __init__(self, trj_file, silent=True):
        self.silent = silent
        self.trj_file = trj_file
        self.data_file = ''
        self.timesteps = []
        self._dump_style = ''
        self._dump_keywords = {} # contains chunk id and keyword. 0 -> 'id'
        self._dump_keywords_r = {} # contains chunk id and keyword. 'id' -> 0
        self.natoms = {}
        self.nheader = 0
        #self.coord = nu.hash() # REMARK: hash seems to be slower than dict in reading trj file
        self.coord = {}
        self.xlo = {}
        self.xhi = {}
        self.ylo = {}
        self.yhi = {}
        self.zlo = {}
        self.zhi = {}
        self.pbc = {}
        self._is_loaded = False
        self._is_dumped = False
        self._nchunk = 0

        # if pickle file already exists, no need to read the file again.
        dir = os.path.dirname(self.trj_file)
        filename = os.path.basename(self.trj_file)
        self.data_file = dir + "." + filename + ".summary.pickle"
        if os.path.exists(self.data_file):
            if not self.silent: nu.warn("Loading information from %s" % self.data_file)
            with open(self.data_file, 'rb') as f:
                timemark = pickle.load(f)
                if timemark == time.ctime(os.path.getmtime(self.trj_file)):
                    if not self.silent: nu.warn("Trajectory file %s is not modified after scan. Proceed to proof-reading.." % self.trj_file)
                    # no modification occured after scan
                    ff = pickle.load(f)
                    self.timesteps = list(set(ff.timesteps))
                    self._dump_style = ff._dump_style
                    self._dump_keywords = ff._dump_keywords
                    self._dump_keywords_r = ff._dump_keywords_r
                    self.natoms = ff.natoms
                    self.nheader = ff.nheader
                    self.xlo = ff.xlo
                    self.xhi = ff.xhi
                    self.ylo = ff.ylo
                    self.yhi = ff.yhi
                    self.zlo = ff.zlo
                    self.zhi = ff.zhi
                    self.pbc = ff.pbc
                    self._is_loaded = True
                    self._nchunk = ff._nchunk

                    return
                else:
                    # file is modified after scan
                    if not self.silent: nu.warn("Trajectory file %s is modified after scan. Reloading the information.." % self.trj_file)

        def scan():
            header = [];
            for line in self._get_line(self.trj_file):
                if "ITEM: ATOMS" in line:
                    header.append(line)
                    self._dump_style = line.strip('\n')
                    dump_keywords = line.replace('ITEM: ATOMS ', '').split()
                    for index, i in enumerate(dump_keywords):
                        self._dump_keywords_r[i] = index
                        self._dump_keywords[index] = i
                    return header
                else:
                    header.append(line)

        scan = scan()
        timestep = int(scan[1])
        self.nheader = len(scan)
        self.natoms[timestep] = int(scan[3])
        self._nchunk = len(scan) + int(scan[3])

    #---------------------------

    def load(self, force=False):
        '''
        Read timestep from trajectory file and stores it to self.timesteps
        '''

        def timesteps_remove_repeat():
            # prevent repeat
            self.timesteps = list(set(self.timesteps))
            self.timesteps.sort()

        if self._is_loaded and not force:
            #nu.warn("Trajectory file already loaded!")
            timesteps_remove_repeat()
            return

        if not self.nheader or not self.natoms:
            nu.warn("Trajectory file %s seems to be empty." % self.trj_file)
            return 0

        dumpatom = self._get_line(self.trj_file)

        while 1:
            try:
                chunk = [next(dumpatom) for i in range(self._nchunk)]
            except StopIteration:
                break;

            t = int(chunk[1])
            self.timesteps.append(t)

            #if not self.silent:
            sys.stdout.write('\rGetting information from LAMMPS trajectory file %s .. Fetched timestep: %d' % (self.trj_file, t))
            sys.stdout.flush()

            self.natoms[t] = int(chunk[3])
            self.xlo[t] = float(chunk[5].split(' ')[0])
            self.xhi[t] = float(chunk[5].split(' ')[1])
            self.ylo[t] = float(chunk[6].split(' ')[0])
            self.yhi[t] = float(chunk[6].split(' ')[1])
            self.zlo[t] = float(chunk[7].split(' ')[0])
            self.zhi[t] = float(chunk[7].split(' ')[1])
            
            self.pbc[t] = [self.xhi[t] - self.xlo[t], self.yhi[t] - self.ylo[t], self.zhi[t] - self.zlo[t]]

        sys.stdout.write('\n')
        sys.stdout.flush()

        # save the information to pickle
        with open(self.data_file, 'wb') as f:
            timemark = time.ctime(os.path.getmtime(self.trj_file))
            pickle.dump(timemark, f)
            pickle.dump(self, f)

        timesteps_remove_repeat()

        self._is_loaded = True

    #---------------------------

    @timer
    def dump(self, requested_ts=[], desc=""):
        if not self._is_loaded:
            nu.warn("LAMMPS trajectory file not loaded. Use load() function first.")
            return

        dumpatom = self._get_line(self.trj_file)

        if not desc: desc = "Dumping trajectories"
        for t in tqdm.tqdm(self.timesteps, ncols=120, desc=desc):
            tinfo = {}
            chunk = [next(dumpatom) for i in range(self._nchunk)]
            if requested_ts:
                if not t in requested_ts:
                    del(chunk)
                    continue

            coords = chunk[9:]
            for line in coords:
                atominfo = {}
                coord = line.split()
                id = int(coord[self._dump_keywords_r['id']])
                type = int(coord[self._dump_keywords_r['type']])
                atominfo['id'] = id
                atominfo['type'] = type
                for index, i in enumerate(coord[2:]):
                    atominfo[self._dump_keywords[index + 2]] = float(i)

                tinfo[id] = atominfo
            self.coord[t] = tinfo

        if len(self.coord) == len(self.timesteps):
            self._is_dumped = True  # True only if timesteps in the trajectory file are fully loaded.

    #---------------------------

    def write(self, target_file):
        with open(target_file, 'w') as f:
            desc = "Writing " + target_file
            timesteps = sorted(self.coord.keys())
            for t in tqdm.tqdm(timesteps, ncols=120, desc=desc):
                # timestep
                f.write("ITEM: TIMESTEP\n")
                f.write("%s\n" % t)

                # number of atoms
                f.write("ITEM: NUMBER OF ATOMS\n")
                f.write("%s\n" % len(self.coord[t]))

                # box bounds    # only supports 3d-periodic trajectory yet: 20160712
                f.write("ITEM: BOX BOUNDS pp pp pp\n")
                f.write(" ".join(str(i) for i in [self.xlo[t], self.xhi[t]]) + "\n")
                f.write(" ".join(str(i) for i in [self.ylo[t], self.yhi[t]]) + "\n")
                f.write(" ".join(str(i) for i in [self.zlo[t], self.zhi[t]]) + "\n")

                # atoms
                f.write(self._dump_style.rstrip() + "\n")
                for atomid in sorted(self.coord[t].keys()):
                    f.write(" ".join([str(self.coord[t][atomid][self._dump_keywords[k]]) for k in self._dump_keywords.keys()]) + "\n")

        return

    #---------------------------

    def split(self, timestep=-1):
        ''' TODO: timestep should be set to i for file
            TODO: if timestep is specified, then the only timestep should be saved.
        '''
        def _chunks(chunk_iterable, n):
            chunk_iterable = iter(chunk_iterable)
            while True:
                yield itertools.chain([next(chunk_iterable)], itertools.islice(chunk_iterable, n-1))

        with open(self.trj_file, 'r') as bigfile:
            for i, lines in enumerate(_chunks(bigfile, self._nchunk)):
                file_split = '{}.{}'.format(self.trj_file, self.timesteps[i])
                with open(file_split, 'w') as f:
                    f.writelines(lines)

        return
