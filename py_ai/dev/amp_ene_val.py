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
from my_arith import rmse
from my_images import Images
import re

def get_title(job, fname, HL, E_conv, ntotal, ndata):
    title = fname.split(".")[0] + "\n"
    hl = '$\\times$'.join(str(x) for x in HL) 
    suptitle = "\n\nAMP Model(HL={}".format(hl) + ", "
    suptitle += "E_rms={}):".format(E_conv) + " "
    if re.search("te", job):
        suptitle += "test/train={}/{}".format(ndata, ntotal-ndata)
    elif re.search("va", job):
        suptitle += "train:validation:test=3:1:1\n"
    return title, suptitle

def exe_train_images(images, HL, E_conv):
    Hidden_Layer=tuple(HL)
    print("Hidden Layer: {}".format(Hidden_Layer))
    print("Energy convergence: {}".format(E_conv))
    calc = Amp(descriptor=Gaussian(), model=NeuralNetwork(hiddenlayers=Hidden_Layer))
    calc.model.lossfunction = LossFunction(convergence={'energy_rmse': E_conv})
    #calc.model.lossfunction = LossFunction(force_coefficient=-0.1)
    calc.train(images=images, overwrite=True)
    return

def run_md(atoms):
    from ase import units
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.md import VelocityVerlet

    traj = ase.io.Trajectory("traj.traj", 'w')

    calc = Amp.load("amp.amp")
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    MaxwellBoltzmannDistribution(atoms, 300. * units.kB)
    traj.write(atoms)
    dyn = VelocityVerlet(atoms, dt=1. * units.fs)
    for step in range(100):
        pot = atoms.get_potential_energy()  # 
        kin = atoms.get_kinetic_energy()
        print("{}: Total Energy={}, POT={}, KIN={}".format(step,pot+kin, pot, kin))
        dyn.run(10)
        traj.write(atoms)

def exe_test_images(job, test_images, amp_pes, title, suptitle,val_id=None):
    calc = Amp.load(amp_pes)
    y=[]
    y_bar=[]
    #return         # this runs
    for mol in test_images:
        y.append(mol.get_potential_energy())
        mol.set_calculator(calc)
        y_bar.append(mol.get_potential_energy())
    return          # this evokes error

    if re.search("te", job):
        draw_dots_two(y, y_bar, title, suptitle)
    elif re.search("va", job):
        err = rmse(y, y_bar)
        print("in job {}-{}: validation error is {}".format(job,val_id,err)) 

def amp_jobs(fdata, job, ndata, HL, E_conv):
    total_images = ase.io.read(fdata, index=':')
    images_c = Images(total_images)
    ### job == training
    if re.search("tr",job):
        if not ndata:
            images = images_c.total_images
            print("Start training using all the data %d" % len(images))
        else:
            images = images_c.get_training_images(ndata)
            print("data training:total sets %d/%d" % (len(images), len(total_images)))
        exe_train_images(images, HL, E_conv)
    ### job == test
    elif re.search("te",job):
        amp_pes = "amp.amp"
        images = images_c.get_test_images(ndata)
        title, suptitle = get_title(job, fdata, HL, E_conv, len(total_images), len(images))
        print("data test:total sets %d/%d" % (len(images), len(total_images)))
        exe_test_images(job, images, amp_pes, title, suptitle)
    ### job == validation
    elif re.search("va",job):
        print("validation test")
        if not ndata:
            ndata = 5
            print("data images are diveded into %d sets" % ndata)
        ### training set scan for valicaiotn
        for i in range(ndata-1):    # last one is kept for test
            ### training
            images, img_valid = images_c.get_val_train_images(ndata, i)
            print("num images: training {} validation {}".format(len(images),len(img_valid)))
            exe_train_images(images, HL, E_conv)
            ### validating
            amp_pes = "amp.amp"
            title, suptitle = get_title(job, fdata, HL, E_conv, len(total_images), len(images))
            #exe_test_images(job, img_valid, amp_pes, title, suptitle, val_id=i)
            ### Alternative Ways
            calc = Amp.load(amp_pes)
            y=[]
            y_bar=[]
            for mol in img_valid:
                y.append(mol.get_potential_energy())
                mol.set_calculator(calc)
                y_bar.append(mol.get_potential_energy())
            '''
            err = rmse(y, y_bar)
            print("in job {}-{}: validation error is {}".format(job,i,err)) 
            '''
            # check divided image sets: plot 2d here
            if False:
                x_draw=[]
                y_draw=[]
                for n, atoms in enumerate(images):
                    pot = atoms.get_potential_energy()
                    x_draw.append(n)
                    y_draw.append(pot)
                mplot_vector_two(x_draw,y_draw, Title="Extracted Training Set %d" % i, Xtitle="serial number", Ytitle="Epot")

    elif re.search('md',job):
        print("ndata is used for start geometry")
        if not ndata:
            atoms = ase.io.read(fdata, index='0')
        else:
            atoms = ase.io.read(fdata, index=ndata)
        run_md(atoms)
    return

def main():
    parser = argparse.ArgumentParser(description='run amp with extxyz ')
    parser.add_argument('fin', help='extxyz input file')
    parser.add_argument('job', default='train', help='job option:"train","test","md","validation"')
    parser.add_argument('-n', '--dlimit', type=int,  help='data range for training and test')
    parser.add_argument('-hl', '--hidden_layer', nargs='*', type=int, default=[8,8,8], help='Hidden Layer of lists of integer')
    parser.add_argument('-el', '--e_convergence', default=0.001, type=float, help='energy convergence limit')
    #group_train = parser.add_mutually_exclusive_group()
    #group_test  = parser.add_mutually_exclusive_group()
    #group_train.add_argument('-l', '--data_limit', type=int, help='the number of data to be trained')
    #group_test.add_argument('-m', '--data_limit2', type=int, help='the start index of data to be tested')
    
    #parser.add_argument('-s', '--sets', default=1, type=int, help='nsets to divide data into training and test set')
    args = parser.parse_args()

    amp_jobs(args.fin, args.job, args.dlimit, args.hidden_layer, args.e_convergence)
    return

if __name__ == '__main__':
    main()

