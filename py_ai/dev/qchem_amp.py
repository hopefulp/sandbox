#!/home/joonho/anaconda3/bin/python

import argparse

import ase.io

from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model.neuralnetwork import NeuralNetwork
from amp.model import LossFunction

def run_amp(fin):
    images = ase.io.read(fin, index='250:')
    print(len(images))

    calc = Amp.load("amp.amp")
    print(images[0])
    for mol in images:
        print(mol.get_potential_energy())
        mol.set_calculator(calc)
        print(mol.get_potential_energy())
    print(images[0])
    
    # For training
    #calc = Amp(descriptor=Gaussian(), model=NeuralNetwork(hiddenlayers=(10, 10, 10)))
    #calc.model.lossfunction = LossFunction(convergence={'energy_rmse': 0.1})
    #calc.model.lossfunction = LossFunction(force_coefficient=-0.1)
    #calc.train(images=images, overwrite=True)

def main():
    parser = argparse.ArgumentParser(description='run amp with extxyz ')
    parser.add_argument('fin', help='extxyz input file')
    args = parser.parse_args()

    run_amp(args.fin)
    return

if __name__ == '__main__':
    main()

