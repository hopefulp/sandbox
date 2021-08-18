#!/home/joonho/anaconda3/bin/python

import argparse
import numpy as np
from ase_build import build_structure
import ase.io
from ase import Atoms, Atom, units
import ase.io
from ase.calculators.lj import LennardJones
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.nptberendsen import NPTBerendsen


def md_run(atoms, count=1000, T=300, trajname="emtmd.traj"):
    traj = ase.io.Trajectory(trajname, 'w')
    atoms.set_calculator(LennardJones(sigma=3.401 ,epsilon=0.001006585 ))  # units: A, eV
    atoms.get_potential_energy()
    traj.write(atoms)
    MaxwellBoltzmannDistribution(atoms, T * units.kB)
    dyn = NPTBerendsen(atoms, timestep=5*units.fs, temperature_K = T, taut=100*units.fs, \
        pressure_au = 0.6889 * units.bar, taup=1000*units.fs, compressibility=980.665 * units.bar) # 
    for step in range(count - 1):
        dyn.run(1)
        traj.write(atoms)
    return 0

def aux_md(inf, job, dyn):
    if inf:
        ### read input geometry
        atoms=ase.io.read(inf)
    else:
        ### use 256 Ar atoms
        atoms = build_structure('bulk', status='gas', natom=256)

    if job == 'traj':
        traj = ase.io.Trajectory(traj_file, "r")
    elif job == 'md':
        md_run(atoms, count=1000, T=80, trajname='Ar-lj.traj')
    return 0    

def main():
    parser = argparse.ArgumentParser(description='analize ASE trajectory')
    parser.add_argument('-inf', '--fin', help='initial geometry for md')
    parser.add_argument('-j', '--job', default='md', help='run ASE md')
    parser.add_argument('-dyn', '--dynamics', default='npt', help='MD scheme')
    args = parser.parse_args()

    aux_md(args.fin, args.job, args.dynamics)
    return

if __name__ == '__main__':
    main()

