#!/home/joonho/anaconda3/bin/python

import argparse

import ase.io
from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model.neuralnetwork import NeuralNetwork
from amp.model import LossFunction

import numpy as np
import matplotlib.pyplot as plt
from my_mplot2d import *

def draw(size, y, h):
    x = range(size)
    fig = plt.figure()
    ax = plt.axes()

    diff =  np.subtract(h,y)
    rmse = np.sqrt((diff**2).mean())
    print(rmse)

    my_font('amp')

    plt.rcParams['figure.figsize'] = (20,10)
    plt.title('AMP model error test-set')
    plt.ylabel('PE(kJ/mol)')
    plt.xlabel('data')
    #plt.scatter(x, y, 'r', x, y_bar, 'b')
    plt.scatter(x, y, marker='^')
    plt.scatter(x, h, marker='o')
    plt.plot(x, diff)
    plt.show()
    

def run_amp(fin):
    images = ase.io.read(fin, index='250:')
    print(len(images))

    calc = Amp.load("amp.amp")
    #print(images[0])
    y=[]
    y_bar=[]
    for mol in images:
        y.append(mol.get_potential_energy())
        #print(mol.get_potential_energy())
        mol.set_calculator(calc)
        y_bar.append(mol.get_potential_energy())
        #print(mol.get_potential_energy())
    #print(images[0])

    draw(len(images), y, y_bar)

    # For training
    #calc = Amp(descriptor=Gaussian(), model=NeuralNetwork(hiddenlayers=(10, 10, 10)))
    #calc.model.lossfunction = LossFunction(convergence={'energy_rmse': 0.1})
    #calc.model.lossfunction = LossFunction(force_coefficient=-0.1)
    #calc.train(images=images, overwrite=True)

    return

def main():
    parser = argparse.ArgumentParser(description='run amp with extxyz ')
    parser.add_argument('fin', help='extxyz input file')
    #parser.add_argument
    args = parser.parse_args()

    run_amp(args.fin)
    return

if __name__ == '__main__':
    main()

