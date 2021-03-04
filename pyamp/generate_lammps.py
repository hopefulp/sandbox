#!/home/joonho/anaconda3/bin/python
import numpy as np
from ase.io import read,write
from amp import Amp
from amp.convert import save_to_prophet
import os
import argparse


def generate_lamphet(amppot, inf, symbols, idx ):
    #Make PROPHET/lammps potentials
    calc = Amp.load(amppot, cores = 1)
    save_to_prophet(calc)


    #atoms = read('starting_configuration.traj')
    atoms = read(inf, index=idx)
    #atoms = images[0]

    # Number of atoms to create
    natoms = len(atoms)
    ntypes = len(set(atoms.get_chemical_symbols()))
    #symbols = ['H', 'O'] #sort them according to atomic number - not alphabetically - to be consistent with masses
    masses = list(set(atoms.get_masses()))
    masses.sort()


    atom_numbers = {}
    for i, symbol in enumerate(symbols):
            atom_numbers[symbol] = i+1

    # System cell
    cell = atoms.get_cell()

    # Generate atom positions
    positions = atoms.get_positions()
        
    # Write LAMMPS data file
    with open('system.data','w') as fdata:
        # First line is a comment line 
        fdata.write('Random atoms - written for EnCodeVentor tutorial\n\n')
        
        #--- Header ---#
        # Specify number of atoms and atom types 
        fdata.write('{} atoms\n'.format(natoms))
        fdata.write('{} atom types\n'.format(ntypes))
        # Specify box dimensions
        fdata.write('{0:16.8e} {1:16.8e} xlo xhi\n'.format(0.0, cell[0,0]))
        fdata.write('{0:16.8e} {1:16.8e} ylo yhi\n'.format(0.0, cell[1,1]))
        fdata.write('{0:16.8e} {1:16.8e} zlo zhi\n'.format(0.0, cell[2,2]))
        fdata.write('{0:16.8e} {1:16.8e} {2:16.8e} xy xz yz\n'.format(cell[0,1], cell[0,2], cell[1,2]))



        # Masses section
        fdata.write(' Masses\n\n')

            # Write each mass
        for i, mass in enumerate(masses):
            fdata.write('{} {}\n'.format(i+1, mass))


        # Atoms section
        fdata.write('\n\n Atoms\n\n')

        # Write each position 
        for i,pos in enumerate(positions):
            fdata.write('{0} {1} {2:16.8e} {3:16.8e} {4:16.8e}\n'.format(i+1, atom_numbers[atoms[i].symbol], pos[0], pos[1], pos[2])) 

def main():
    parser = argparse.ArgumentParser(description='Process for Genetic Algorithm & Artificial Neural Network')
    parser.add_argument('-p', '--amppot', default='amp.amp', help='amp potential')
    parser.add_argument('-f', '--inf', default='OUTCAR', help='input file')
    parser.add_argument('-s', '--symbols', nargs='+', default=['H','O'], help='atom list')
    parser.add_argument('-i', '--index', default=0, help='get the initial structure from input file')
    args = parser.parse_args()

    generate_lamphet(args.amppot, args.inf, args.symbols, args.index)


if __name__ == '__main__':
    main()
