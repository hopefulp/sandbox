#!/home/joonho/anaconda3/bin/python
'''
    2020.08.25 GA(genetic algorithm) was encoded by job='trga', 'tega', if 'ga' in job, turn on Lga
    2020.08.25 Ltest-force is deprecated. if there is force training, make a force test
    2020.11.13 find amp-pot in the directory even though -p None
'''
import argparse
from amp import Amp
import ase
from ase import units
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md import VelocityVerlet

def amp_md(atoms, nstep, dt, amp_pot):

    traj = ase.io.Trajectory("traj.traj", 'w')
    calc = Amp.load(amp_pot)

    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    MaxwellBoltzmannDistribution(atoms, 300 * units.kB)
    traj.write(atoms)
    dyn = VelocityVerlet(atoms, dt=dt * units.fs)
    f = open("aimd.ene", "w")
    f.write(f"{'time':^5s}{'Etot':^15s}{'Epot':^15s}{'Ekin':^10s}\n")
    print(f"   {'time':^5s}{'Etot':^15s}{'Epot':^15s}{'Ekin':^10s}")
    dump_freq = 50
    nstep = nstep/dump_freq
    for step in range(int(nstep)):
        pot = atoms.get_potential_energy()  # 
        kin = atoms.get_kinetic_energy()
        tot = pot + kin
        f.write(f"{step:5d}{tot:15.4f}{pot:15.4f}{kin:10.4f}\n")
        print(f"{step:5d}{tot:15.4f}{pot:15.4f}{kin:10.4f}")
        dyn.run(dump_freq)
        traj.write(atoms)                   # write kinetic energy, but pot is not seen in ase
    f.close()        

def main():
    parser = argparse.ArgumentParser(description='run amp with extxyz, OUTCAR: validation is removed ', prefix_chars='-+/')
    parser.add_argument('-p', '--pot', default='amp.amp', help="input amp potential")
    parser.add_argument('-fin', '--infile', help='input trajectory file which can be read by ase')
    parser.add_argument('-i','--index', default=0, type=int, help='select start configuration from input file')
    parser.add_argument('-ns','--nstep', default=1000, type=int, help='number of step with dt')
    parser.add_argument('-dt','--dt', default=1.0, type=float, help='time interval in fs')

    args = parser.parse_args()

    atoms = ase.io.read(args.infile, index=args.index)      # index = string
    amp_md(atoms, args.nstep, args.dt, args.pot)

if __name__ == '__main__':
    main()

