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
from my_arith import divide_int
"""
fig = plt.figure(figsize=(15,10))
ax = plt.axes()
text_pos_x = 55
text_pos_y = 10

def draw(y, h, title, suptitle):

    nlen = len(y)
    h_conv = np.array(h) * ev2kj
    y_conv = np.array(y) * ev2kj
    diff =  np.subtract(h_conv,y_conv)
    rmse = np.sqrt((diff**2).mean())
    max_res = abs(max(diff, key=abs))
    #max_res = max(diff, key=abs)
    #print("{:10.3f} {:10.3f}".format(rmse,max_res))
    ### input text inside figure
    text="E_rms(test) = {:7.3f}\nE_maxres = {:7.3f}".format(rmse, max_res)
    plt.text(text_pos_x, text_pos_y,text, fontsize=20)


    ones = np.zeros((len(y_conv)))
    #my_font('amp')
    mpl.rcParams.update({'font.size':22})
    plt.title(title, fontsize=40)
#    plt.suptitle(suptitle, x=0.5, y=0.92, va='top', fontsize=18)
    plt.suptitle(suptitle, fontsize=18)
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
"""
   
 
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

def divide_train_valid_images(images_train_all, idx_imagelist, n):
    '''
    images_train_all:: exclude the test images
    idx_imagelist:: divided index is up to last sets for validation
    '''
    ndata = len(idx_imagelist) + 1

    if n == 0:
        v_start = 0
        v_stop  = idx_imagelist[0]
    else:
        v_start = idx_imagelist[n-1]
        v_stop  = idx_imagelist[n]
    valid_images=images_train_all[v_start:v_stop]
    if n == 0:
        train_images=images_train_all[v_stop:idx_imagelist[3]]
    else:
        train_images=images_train_all[0:v_start]
        if n != ndata-2:
            train_images.extend(images_train_all[v_stop:idx_imagelist[3]])
    return train_images, valid_images

def exe_test_images(images_4test, amp_pes, title, suptitle):
    calc = Amp.load(amp_pes)
    y=[]
    y_bar=[]
    for mol in images_4test:
        y.append(mol.get_potential_energy())
        mol.set_calculator(calc)
        y_bar.append(mol.get_potential_energy())
    draw_dots_two(y, y_bar, title, suptitle)

def amp_jobs(fdata, job, ndata, HL, E_conv):
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
        exe_train_images(images, HL, E_conv)

    elif job == "test":
        amp_pes = "amp.amp"

        ### get model and make title
        title = fdata.split(".")[0] + "\n"
        hl = '$\\times$'.join(str(x) for x in HL) 
        suptitle = "\n\nAMP Model(HL={}".format(hl) + ", "
        suptitle += "E_rms={}):".format(E_conv) + " "
        suptitle += "test/train={}/{}".format(len(total_images)-ndata, ndata)

        exe_test_images(total_images[ndata:], amp_pes, title, suptitle)

    elif job == "valid-test":
        print("validation test")
        if not ndata:
            ndata = 5
        parts = divide_int(len(total_images), ndata)
        print(parts, "in total {}".format(len(total_images)))
        for i in range(ndata-1):
            img_train, img_valid = divide_train_valid_images(total_images[0:parts[-1]], parts, i)
            
            exe_train_images(img_train, HL, E_conv)
            amp_pes = "amp.amp"
            ### get model and make title
            title = fdata.split(".")[0] + "\n"
            hl = '$\\times$'.join(str(x) for x in HL)
            suptitle = "AMP Model(HL={}".format(hl) + ", "
            suptitle += "E_rms={}):".format(E_conv) + " "
            suptitle += "train:validation:test=1:1:1\n"
            exe_test_images(total_images[ndata:], amp_pes, title, suptitle)

            # check divided image sets: plot 2d here
            if False:
                x_draw=[]
                y_draw=[]
                for n, atoms in enumerate(img_train):
                    pot = atoms.get_potential_energy()
                    x_draw.append(n)
                    y_draw.append(pot)
                mplot_vector_two(x_draw,y_draw, Title="Extracted Training Set %d" % i, Xtitle="serial number", Ytitle="Epot")

    elif job == 'md':
        if not ndata:
            atoms = ase.io.read(fdata, index='0')
        else:
            atoms = ase.io.read(fdata, index=ndata)
        run_md(atoms)
    return

def main():
    parser = argparse.ArgumentParser(description='run amp with extxyz ')
    parser.add_argument('fin', help='extxyz input file')
    parser.add_argument('job', default='train', choices=['train','test','md','valid-test'], help='job option: "train", "test", "md"')
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

