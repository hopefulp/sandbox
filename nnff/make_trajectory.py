#!/home/joonho/anaconda3/bin/python
from ase.io import read,write
from ase.visualize import view
from ase.calculators.singlepoint import SinglePointCalculator
from ase import units
import numpy as np
import argparse

#Get 1ps in ASE units
ps = 1000*units.fs

def make_ase_traj(atom_symbols):

    #Atoms present in simulation (needed to convert lammps atom types to the correct atoms in ASE)
    #atom_symbols = ['H', 'O']


    #read what is possible from lammps dump file
    # ase 3.19.1: use format=lammps-dump-text
    # ase 3.19.3 <= : makes an error
    try:
        images = read('dump.lammps', format='lammps-dump-text', index = ':')
        Nsteps = len(images) #we will use this to know how many energy lines to read
    except IndexError:
        print('Atoms lost in last steps... Trajectory will not contain last two structures.')
        images = read('dump.lammps', format='lammps-dump-text', order = False, index = ':')
        images = images[0:-2]
        Nsteps = len(images)



    #Read energies
    logfile =  open('log.lammps', 'r')
    lines = logfile.readlines()

    #Get spacing between each dump. This is probably a pretty stupid way to get it, but it works... 
    for line in lines[0:1000]:
        if 'dump ' in line:
            spacing = int(line.split('custom ')[1].split(' dump.lammps')[0]) #the dump spacing is sandwiched between the words "format" and "dump.lammps"
            break


    epot = []
    ekin = []

    dump_file = open('dump.lammps', 'r')
    first_timestep = int(dump_file.readlines()[1])
    index = lines.index('Step PotEng KinEng \n')+1+first_timestep
    energy_lines = lines[index:index+Nsteps*spacing:spacing]
    for line in energy_lines:

        line = (line.strip().split(' '))
        line = list(filter(('').__ne__, line))
        
        try:
            epot.append(float(line[1]))
        except ValueError: 
            break
        ekin.append(float(line[2]))


    #Read forces and correct atomic symbols and velocites
    for i, atoms in enumerate(images):
        atomic_numbers = atoms.get_atomic_numbers()
        forces = atoms.get_forces()
        try:
            energy = epot[i]
        except IndexError:
            energy = epot[i-1]
        velocities = atoms.get_velocities()
        velocities = velocities/ps #convert velocities to ASE units

        for j, atom in enumerate(atoms):
            atom.symbol = atom_symbols[atomic_numbers[j]-1]

        atoms.set_pbc([True,True,True])
        calc = SinglePointCalculator(atoms, energy = energy, forces = forces)
        atoms.set_calculator(calc)

        #Set velocities
        atoms.set_velocities(velocities)

    write('MD_trajectory.traj', images)

def main():
    parser = argparse.ArgumentParser(description='Process for Genetic Algorithm & Artificial Neural Network')
    parser.add_argument('-s', '--symbols', nargs='+', default=['H','O'], help='atom list')
    args = parser.parse_args()

    make_ase_traj(args.symbols)


if __name__ == '__main__':
    main()

