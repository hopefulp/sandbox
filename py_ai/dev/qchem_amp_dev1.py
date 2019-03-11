#!/home/joonho/anaconda3/bin/python

import argparse

import ase.io

from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model.neuralnetwork import NeuralNetwork
from amp.model import LossFunction

def run_amp(fin):
    images = ase.io.read(fin, index=':')

    print(len(images))

    calc = Amp(descriptor=Gaussian(), model=NeuralNetwork(hiddenlayers=(10, 10, 10)))
    calc.model.lossfunction = LossFunction(convergence={'energy_rmse': 0.001})
    calc.model.lossfunction = LossFunction(force_coefficient=-0.1)
    calc.train(images=images, overwrite=True)

def main():
    parser = argparse.ArgumentParser(description='run amp with extxyz and low ')
    parser.add_argument('fin', help='extxyz input file')
    args = parser.parse_args()

    run_amp(args.fin)
    return

if __name__ == '__main__':
    main()

