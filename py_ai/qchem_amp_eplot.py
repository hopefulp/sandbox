#!/home/joonho/anaconda3/bin/python

import argparse

import ase.io
from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model.neuralnetwork import NeuralNetwork
from amp.model import LossFunction

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from my_mplot2d import *
from my_chem import *

fig = plt.figure(figsize=(15,10))
ax = plt.axes()

def draw(y, h):

    nlen = len(y)
    h_conv = np.array(h) * ev2kj
    y_conv = np.array(y) * ev2kj
    diff =  np.subtract(h_conv,y_conv) 
    rmse = np.sqrt((diff**2).mean())
    max_res = abs(max(diff, key=abs))
    #max_res = max(diff, key=abs)
    print("{:10.3f} {:10.3f}".format(rmse,max_res))
    ones = np.zeros((len(y_conv)))

    #my_font('amp')
    mpl.rcParams.update({'font.size':22})
    plt.title('AMP model error test-set')
    plt.ylabel('PE(kJ/mol)', fontsize=20)
    plt.xlabel('data', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=15)
    #plt.scatter(x, y, 'r', x, y_bar, 'b')
    p1  = plt.scatter(range(nlen), y_conv, c='r', marker='o', label='true value')
    p2  = plt.scatter(range(nlen), h_conv, c='b', marker='^', label='hypothesis')
    p3, = plt.plot(range(nlen), diff, label='difference')
    plt.plot(range(nlen), ones)
    plt.legend([p1,p2,p3],['true value', 'hypothesis', 'difference'])
    plt.show()

    return
    
def train_images(images):
    calc = Amp(descriptor=Gaussian(), model=NeuralNetwork(hiddenlayers=(10, 10, 10)))
    calc.model.lossfunction = LossFunction(convergence={'energy_rmse': 0.001})
    calc.model.lossfunction = LossFunction(force_coefficient=-0.1)
    calc.train(images=images, overwrite=True)


def run_amp(fdata, job, ndata):

    total_images = ase.io.read(fdata, index=':')
    if job == "train":
        if not ndata:
            #images = ase.io.read(fdata, index=':')
            images = total_images
            print("Start training using all the data %d" % len(total_images))
        else:
            #images = ase.io.read(fdata, index=':ndata')
            images = total_images[:ndata]
            print("number of traing data %d/%d" % (len(images), len(total_images)))
        train_images(images)

    elif job == "test":
        #images = ase.io.read(fdata, index='ndata:')
        images = total_images[ndata:]

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

        draw(y, y_bar)
    return

def main():
    parser = argparse.ArgumentParser(description='run amp with extxyz ')
    parser.add_argument('fin', help='extxyz input file')
    parser.add_argument('job', default='train', choices=['train', 'test'], help='job option: "train" vs "test"')
    parser.add_argument('-n', '--dlimit', type=int,  help='data range for training and test')
    group_train = parser.add_mutually_exclusive_group()
    group_test  = parser.add_mutually_exclusive_group()
    #group_train.add_argument('-l', '--data_limit', type=int, help='the number of data to be trained')
    #group_test.add_argument('-m', '--data_limit2', type=int, help='the start index of data to be tested')
    
    #parser.add_argument('-s', '--sets', default=1, type=int, help='nsets to divide data into training and test set')
    args = parser.parse_args()

    run_amp(args.fin, args.job, args.dlimit)
    return

if __name__ == '__main__':
    main()

