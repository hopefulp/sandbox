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
    
def train_all(images):
    calc = Amp(descriptor=Gaussian(), model=NeuralNetwork(hiddenlayers=(10, 10, 10)))
    calc.model.lossfunction = LossFunction(convergence={'energy_rmse': 0.001})
    calc.model.lossfunction = LossFunction(force_coefficient=-0.1)
    calc.train(images=images, overwrite=True)


def run_amp(fdata, job, ndata):


    """
    if nsets == 1:
        train_all(images)
    else:
        job = 'test'
        ndata = len(images)
        data_sets = divide_dataset(ndata, nsets)

        images = ase.io.read(fin, index='250:')
    """
    print(ndata)

    if job == "all":
        images = ase.io.read(fdate, index=':')
        train_all(images)
    elif job == "train":
        images = ase.io.read(fdate, index=':ndata')

    ndata = len(images)

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

    return

def main():
    parser = argparse.ArgumentParser(description='run amp with extxyz ')
    parser.add_argument('fin', help='extxyz input file')
    parser.add_argument('job', default='all', choices=['all', 'train', 'test'], help='job option: \n all: train all the data \n train: train part of data \n test: test the trained model')
    parser.add_argument('-n', '--dlimit', type=int, help='data range for training and test')
    group_all   = parser.add_mutually_exclusive_group()
    group_train = parser.add_mutually_exclusive_group()
    group_test  = parser.add_mutually_exclusive_group()
    #group_train.add_argument('-l', '--data_limit', type=int, help='the number of data to be trained')
    #group_test.add_argument('-m', '--data_limit2', type=int, help='the start index of data to be tested')
    
    #parser.add_argument('-s', '--sets', default=1, type=int, help='nsets to divide data into training and test set')
    args = parser.parse_args()

    run_amp(args.fin, args.job, args.nlimit)
    return

if __name__ == '__main__':
    main()

