#!/home/joonho/anaconda3/bin/python

from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units

from ase.calculators.emt import EMT
size = 5

import ase.io
from ase.io.trajectory import Trajectory

atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                         symbol='Pd',
                         size=(size, size, size),
                         pbc=True)

traj = ase.io.Trajectory('traj_1.traj','w')

atoms.set_calculator(EMT())

MaxwellBoltzmannDistribution(atoms, 300 * units.kB)

traj.write(atoms)

dyn = VelocityVerlet(atoms, dt=1 * units.fs)

for step in range(10):
    pot = atoms.get_potential_energy()
    kin = atoms.get_kinetic_energy()
    with open('traj_1.txt','w') as f:
        f.write("{}: Total Energy={}, POT={}, KIN={}\n".format(step, pot+kin, pot, kin))
    dyn.run(5)
    ase.io.write('traj_1.xyz',ase.io.read('traj_1.traj'), append=True)
    traj.write(atoms)


