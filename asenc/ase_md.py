#!/home/joonho/anaconda3/bin/python

import argparse
import numpy as np

import ase.io
from ase import Atoms, Atom, units
import ase.io
from ase.calculators.emt import EMT
from ase.build import fcc110
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md import VelocityVerlet
from ase.constraints import FixAtoms


def md_run(inf, count=1000, filename="emtmd.traj"):
    atoms=ase.io.read(inf)
    traj = ase.io.Trajectory(filename, 'w')
    atoms.set_calculator(EMT())
    atoms.get_potential_energy()
    traj.write(atoms)
    MaxwellBoltzmannDistribution(atoms, 300. * units.kB)
    dyn = VelocityVerlet(atoms, dt=1. * units.fs)
    for step in range(count - 1):
        dyn.run(1)
        traj.write(atoms)
    return 0

def anal_traj(traj_file, job):
    if job == 'traj':
        traj = ase.io.Trajectory(traj_file, "r")
    elif job == 'md':
        md_run(traj_file)
    return 0    

def main():
    parser = argparse.ArgumentParser(description='analize ASE trajectory')
    parser.add_argument('fin', help='ASE traj file')
    parser.add_argument('job', default='traj', choices=['traj','md','ampmd'], help='job option: "traj"')
    args = parser.parse_args()

    anal_traj(args.fin, args.job)
    return

if __name__ == '__main__':
    main()

