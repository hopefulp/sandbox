#!/usr/bin/python

import argparse
import re
import os
import numpy as np

from ase import Atoms, Atom
import ase.io

from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model.neuralnetwork import NeuralNetwork

def amp_run(data_file, fit_type):

    ### Data obtain
    data = np.load(d_file)
    energies = data['ene']
    coordinates = data['coord']

    traj = ase.io.Tr

    return

def main():
    parser = argparse.ArgumentParser(description='LOAD data file (npz), RUN AMP')
    parser.add_argument('QM_data', help='Load data file (npz)')
    parser.add_argument('-t', '--fit_type', default='e', help='data type: energy alone, energy & force, energy, force & hessian')
    args = parser.parse_args()

    amp_run(args.QM_data, args.fit_type) 

if __name__ == '__main__':
    main()
