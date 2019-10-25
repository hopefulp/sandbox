#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys
import os
import math
import itertools
import copy
import numpy as np

import tqdm
import memory_footprint
import lammpstrj as lt
import datatools as dt

version = "160715"
usage = """Usage: %s dat_file trj_file mode
Creates a series of LAMMPS trajectory files with 
1) vz suppressed trajectory (vz = 0)
2) vx&vy suppressed trajectory (vx = vy = 0)
3) Vcm suppressed & vz suppressed trajectory 
4) Vcm suppressed & vx&vy suppressed trajectory 

""" % sys.argv[0]

def myhandler(obj):
    assert isinstance(obj, lt.lammpstrj)
    yield obj.trj_file
    yield obj.data_file
    yield obj.timesteps
    yield obj._dump_style
    yield obj._dump_keywords
    yield obj._dump_keywords_r
    yield obj.natoms
    yield obj.nheader
    yield obj.coord
    yield obj.xlo
    yield obj.xhi
    yield obj.ylo
    yield obj.yhi
    yield obj.zlo
    yield obj.zhi
    yield obj.pbc
    yield obj._is_loaded
    yield obj._is_dumped

def chunk(it, n):
    it = iter(it)
    return iter(lambda: list(itertools.islice(it, n)), [])

if len(sys.argv) < 1:
    print(usage)
    sys.exit()

# options
dat_file = sys.argv[1]
trj_file = sys.argv[2]
try:
    mode = int(sys.argv[3])
except:
    mode = 1

# output files
'''
if mode == 1:
    out_trj_file = trj_file.split('.lammps')[0] + '.vz0.lammps'
elif mode == 2:
    out_trj_file = trj_file.split('.lammps')[0] + '.vxy0.lammps'
elif mode == 3:
    out_trj_file = trj_file.split('.lammps')[0] + '.Vcm-vz0.lammps'
elif mode == 4:
    out_trj_file = trj_file.split('.lammps')[0] + '.Vcm-vxy0.lammps'
else:
    nu.die("Wrong mode specified.")
'''
out_trj_file = trj_file.split('.lammps')[0] + '.mod.lammps'

print("* Requested mode: %d" % mode)
print("* Modified trajectory will be saved at %s" % out_trj_file)

# read masses and connections from data file
print("Reading atom masses and connections from data file..")
masses = dt.get_mass(dat_file)
conn = dt.get_connection(dat_file)
atom_types = dt.get_atom_types(dat_file)


# dump coords & velocities from the trajectory
mytrj = lt.lammpstrj(trj_file)
mytrj.load()
timesteps = sorted(mytrj.timesteps)


# estimate memory usage and slice dump file 
mytrj.dump(requested_ts=[timesteps[0]], desc="Scanning")   # contains only one shot
avail_memory = memory_footprint.get_avail_memory()['free'] * 1024 * 0.8    # in Byte, 0.8 is the factor of the reserved memory for the script
sizeof_oneshot = memory_footprint.getsizeof(mytrj, {lt.lammpstrj: myhandler})   # in Byte
sizeof_dump = sizeof_oneshot * len(timesteps)
n_pass = math.ceil(sizeof_dump / avail_memory)


# print memory stats
print("* Available memory: %s" % memory_footprint.sizeof_fmt(avail_memory))
print("* Estimated memory usage for one timestep: %s" % memory_footprint.sizeof_fmt(sizeof_oneshot))
print("* Estimated memory usage for whole trajectory: %s" % memory_footprint.sizeof_fmt(sizeof_dump))
print("* # of passes for dump: %d" % n_pass)

n_timesteps_per_pass = math.ceil(float(len(timesteps)) / float(n_pass))
pass_timesteps = list(chunk(timesteps, n_timesteps_per_pass))
print("* # of timesteps for one pass: %d" % n_timesteps_per_pass)
print("* Estimated memory usage for one pass: %s" % memory_footprint.sizeof_fmt(sizeof_oneshot * n_timesteps_per_pass))


# sweep
pass_trj_files = []
for index, pass_t in enumerate(pass_timesteps):
    pass_t = sorted(pass_t)
    pass_num = index + 1
    print("\n=== Pass %d / %d ===" % (pass_num, len(pass_timesteps)))

    # dump trajectory
    mytrj = lt.lammpstrj(trj_file)
    mytrj.dump(requested_ts=pass_t)

    # suppress
    if mode == 1:
        for t in tqdm.tqdm(pass_t, ncols=120, desc='Suppressing vz'):
            for atom in mytrj.coord[t].keys():
                mytrj.coord[t][atom]['vz'] = 0.0

        pass_trj_file = out_trj_file + ".%d" % pass_num

    elif mode == 2:
        for t in tqdm.tqdm(pass_t, ncols=120, desc='Suppressing vx & vy'):
            for atom in mytrj.coord[t].keys():
                mytrj.coord[t][atom]['vx'] = 0.0
                mytrj.coord[t][atom]['vy'] = 0.0
        pass_trj_file = out_trj_file + ".%d" % pass_num

    elif mode == 3:
        for t in tqdm.tqdm(pass_t, ncols=120, desc='Removing Vcm & vz->0'):
            for molecule in conn:
                # calculate Vcm = (m1v1 + m2v2 + ...) / M
                mivi = np.array([0.0, 0.0, 0.0]); M = 0.0
                for i in molecule:
                    vi = np.array([mytrj.coord[t][i]['vx'], mytrj.coord[t][i]['vy'], mytrj.coord[t][i]['vz']])
                    mivi += masses[atom_types[i]] * vi
                    M += masses[atom_types[i]]
                v_cm = mivi / M

                # substract & suppress
                for i in molecule:
                    mytrj.coord[t][i]['vx'] -= vi[0]
                    mytrj.coord[t][i]['vy'] -= vi[1]
                    mytrj.coord[t][i]['vz'] = 0
        pass_trj_file = out_trj_file + ".%d" % pass_num

    elif mode == 4:
        for t in tqdm.tqdm(pass_t, ncols=120, desc='Removing Vcm & vx&vy->0'):
            for molecule in conn:
                # calculate Vcm = (m1v1 + m2v2 + ...) / M
                mivi = np.array([0.0, 0.0, 0.0]); M = 0.0
                for i in molecule:
                    vi = np.array([mytrj.coord[t][i]['vx'], mytrj.coord[t][i]['vy'], mytrj.coord[t][i]['vz']])
                    mivi += masses[atom_types[i]] * vi
                    M += masses[atom_types[i]]
                v_cm = mivi / M

                # substract & suppress
                for i in molecule:
                    mytrj.coord[t][i]['vx'] = 0
                    mytrj.coord[t][i]['vy'] = 0
                    mytrj.coord[t][i]['vz'] -= vi[2]
        pass_trj_file = out_trj_file + ".%d" % pass_num

    else:
        nu.die("Wrong mode number specified.")

    # write trajectory
    mytrj.write(pass_trj_file)
    pass_trj_files.append(pass_trj_file)


# collect trajectories
pass_trj_files2 = " ".join(str(i) for i in pass_trj_files)
print("Concatenating..")
cmd_concat = "cat " + pass_trj_files2 + " >> " + out_trj_file
os.system(cmd_concat)
print("- Multipass files are concatenated to %s" % out_trj_file)


# remove temporary trajectories (*.1, *.2, ...)
print("Removing temp trajectories..")
for i in pass_trj_files:
    print("- Removing %s" % i)
    os.remove(i)


# Done.
print("Done.")
