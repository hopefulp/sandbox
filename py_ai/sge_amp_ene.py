#!/gpfs/home/joonho/anaconda3/bin/python

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
import my_chem
from my_arith import rmse
from my_images import Images
import re
import sys

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
    #calc.model.lossfunction = LossFunction(convergence={'energy_rmse': E_conv},force_coefficient=-0.1)
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
    MaxwellBoltzmannDistribution(atoms, 10. * units.kB)
    traj.write(atoms)
    dyn = VelocityVerlet(atoms, dt=1. * units.fs)
    f = open("md.ene", "w")
    f.write("{:^5s}{:^10s}{:^10s}{:^10s}\n".format("time","Etot","Epot","Ekin"))
    for step in range(100):
        pot = atoms.get_potential_energy()  # 
        kin = atoms.get_kinetic_energy()
        tot = pot + kin
        f.write("{:5d}{:10.5f}{:10.5f}{:10.5f}\n".format(step, tot, pot, kin))
        print("{}: Total Energy={}, POT={}, KIN={}".format(step, tot, pot, kin))
        dyn.run(10)
        traj.write(atoms)
    f.close()        

def exe_test_images(job, test_images, amp_pes, title, suptitle,Lgraph,val_id=None):
    calc = Amp.load(amp_pes)
    y=[]
    y_bar=[]
    #return         # this runs
    for mol in test_images:
        y.append(mol.get_potential_energy())
        mol.set_calculator(calc)
        y_bar.append(mol.get_potential_energy())

    if Lgraph:
        err = draw_dots_two(y, y_bar, title, suptitle)
    else:
        err = rmse(y, y_bar)*my_chem.ev2kj
    return err

def amp_jobs(fdata, job, nsets, HL, E_conv, Lgraph,ival_set):
    total_images = ase.io.read(fdata, index=':')
    images_sets = Images(total_images, nsets)
    if re.search("pr", job):
        y=[]
        for mol in total_images:
            y.append(mol.get_potential_energy())
        mplot_nvector([],y,fdata.split(".")[0],'sample','E(eV)')
    ### job == training
    elif re.search("tr",job):
        images = images_sets.get_training_images()
        print("data training:total sets %d/%d" % (len(images), len(total_images)))
        exe_train_images(images, HL, E_conv)
    ### job == training & test - test can be done at once by commenting one line below
        amp_pes = "amp.amp"
        images = images_sets.get_test_images()
        title, suptitle = get_title(job, fdata, HL, E_conv, len(total_images), len(images))
        print("data test:total sets %d/%d" % (len(images), len(total_images)))
        exe_test_images(job, images, amp_pes, title, suptitle,Lgraph)
    ### only test
    elif re.search("te",job):
        amp_pes = "amp.amp"
        images = images_sets.get_test_images()
        title, suptitle = get_title(job, fdata, HL, E_conv, len(total_images), len(images))
        print("data test:total sets %d/%d" % (len(images), len(total_images)))
        exe_test_images(job, images, amp_pes, title, suptitle,Lgraph)
    ### job == validation
    elif re.search("va",job):
        print("validation test")
        print("data images are diveded into %d sets" % nsets)
        ### training set scan for valicaiotn
        #for i in [0,1,2,3]:   #range(nsets-1):    # last one [4] is kept for test
        # ival_set should be lower than nsets-1
        
        if ival_set is None:
            print("index for validation set is reguired with '-i num' between 0 ~ {}".format(nsets-2))
            sys.exit(0)
        else:
            if ival_set >= nsets-1:
                print("validation set index should be lower than {}".format(nsets-1))
                print("refer to py_ai_ini.py -j amp")
                sys.exit(3)
        fname = fdata.split(".")[0]
        hl = ''.join(str(x) for x in HL)
        fname += hl + str(E_conv) + ".val"
        #for i in range(nsets-1): # last one is kept for test, this is not working at the moment
        for i in [ival_set]:   
            ### training
            images, img_valid = images_sets.get_val_train_images(i)
            print("num images: training {} validation {}".format(len(images),len(img_valid)))
            exe_train_images(images, HL, E_conv)
            ### validating
            amp_pes = "amp.amp"
            title, suptitle = get_title(job, fdata, HL, E_conv, len(total_images), len(images))
            rmserr = exe_test_images(job, img_valid, amp_pes, title, suptitle, Lgraph,val_id=i)
            with open(fname, "a") as f:
                f.write("{}: {:5.3f}\n".format(ival_set,rmserr))
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
        #print("nsets is used for start geometry")
        if not nsets:
            atoms = ase.io.read(fdata, index='0')
        else:
            atoms = ase.io.read(fdata, index=nsets)
        run_md(atoms)
    return

def main():
    parser = argparse.ArgumentParser(description='run amp with extxyz ', prefix_chars='-+/')
    parser.add_argument('fin', help='extxyz input file')
    parser.add_argument('job', default='train', help='job option:"train","test","md","validation","profile"')
    parser.add_argument('-n','--nsets',default=5,type=int,help='num of sets:1 train all sets, otherwise, last set is for test')
    parser.add_argument('-hl', '--hidden_layer', nargs='*', type=int, default=[8,8,8], help='Hidden Layer of lists of integer')
    parser.add_argument('-el', '--e_convergence', default=0.001, type=float, help='energy convergence limit')
    parser.add_argument('-g', action="store_false", help='if val default is False, otherwise True')
    parser.add_argument('+g', action="store_true", help='if val default is False, otherwise True')
    group_valid = parser.add_argument_group()
    group_valid.add_argument('-i', '--index_val_set', type=int, help='validation set index')
    args = parser.parse_args()

    if re.search("tr", args.job):
        args.g = True
    amp_jobs(args.fin, args.job, args.nsets, args.hidden_layer, args.e_convergence,args.g,args.index_val_set)
    return

if __name__ == '__main__':
    main()

