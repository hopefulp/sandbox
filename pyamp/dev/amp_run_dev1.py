#!/home/joonho/anaconda3/bin/python

import argparse

import ase.io
from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model.neuralnetwork import NeuralNetwork
from amp.model import LossFunction

import numpy as np
#from myplot2D import *                 ### FOR MLET 
import my_chem
from my_arith import rmse
from my_images import Images
from common import yes_or_no
from amp_plot import get_title         ### FOR MLET
import re
import sys
import os
import socket

'''
# move to module amp_plot
def get_title(job, fname, HL, E_conv,f_conv, ntotal, ndata):
    title = fname.split(".")[0] + "\n"
    hl = '$\\times$'.join(str(x) for x in HL) 
    suptitle = "\n\nAMP Model(HL={}".format(hl) + ", "
    suptitle += "E_rms={} ".format(E_conv) + " "
    suptitle += "F_err={}):".format(f_conv) + " "
    if re.search("te", job):
        suptitle += "test/train={}/{}".format(ndata, ntotal-ndata)
    elif re.search("va", job):
        suptitle += "train:validation:test=3:1:1\n"
    return title, suptitle
'''
def calc_train_images(images, HL, E_conv, f_conv, ncore):
    Hidden_Layer=tuple(HL)
    print("Hidden Layer: {}".format(Hidden_Layer))
    print("Energy convergence: {}".format(E_conv))
    calc = Amp(descriptor=Gaussian(), model=NeuralNetwork(hiddenlayers=Hidden_Layer), cores=ncore)
    if f_conv <= 0.0:
        calc.model.lossfunction = LossFunction(convergence={'energy_rmse': E_conv})
    else:
        calc.model.lossfunction = LossFunction(convergence={'energy_rmse': E_conv,'force_rmse':f_conv})
        
    #calc.model.lossfunction = LossFunction(force_coefficient=-0.1)
    calc.train(images=images, overwrite=True)
    return

def run_md(atoms):
    from ase import units
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.md import VelocityVerlet

    traj = ase.io.Trajectory("traj.traj", 'w')

    try:
        calc = Amp.load("amp.amp")
    except FileNotFoundError:
        try:
            calc = Amp.load("amp-untrained-parameters.amp") 
        except FileNotFoundError:
            print("Error: amp-pes.amp file does not exist, input amp-pot file by -p")
            sys.exit(1)

    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    MaxwellBoltzmannDistribution(atoms, 300 * units.kB)
    traj.write(atoms)
    dyn = VelocityVerlet(atoms, dt=1. * units.fs)
    f = open("md.ene", "w")
    f.write("{:^5s}{:^10s}{:^10s}{:^10s}\n".format("time","Etot","Epot","Ekin"))
    for step in range(3):
        pot = atoms.get_potential_energy()  # 
        kin = atoms.get_kinetic_energy()
        tot = pot + kin
        f.write("{:5d}{:10.5f}{:10.5f}{:10.5f}\n".format(step, tot, pot, kin))
        print("{}: Total Energy={}, POT={}, KIN={}".format(step, tot, pot, kin))
        dyn.run(10)
        traj.write(atoms)
    f.close()        

def calc_test_images(job, test_images, amp_pes, title, suptitle,ncore, Lgraph,val_id=None,nmol=1,Ltwinx=None):
    try:
        calc = Amp.load(amp_pes)
    except FileNotFoundError:
        try:
            calc = Amp.load("amp-untrained-parameters.amp")
        except FileNotFoundError:
            print("Error: amp-pes.amp file does not exist, input amp-pot file by -p")
            sys.exit(1)
    escale = 1                  # my_chem.ev2kj
    y=[]
    y_bar=[]
    #return         # this runs
    for mol in test_images:
        y.append(mol.get_potential_energy()/nmol)
        mol.set_calculator(calc)
        y_bar.append(mol.get_potential_energy()/nmol)
    
    with open("amp_test.txt",'w') as f:
        f.write("\ttrue\t\thypothesis\n")
        for y1, ybar1 in zip(y, y_bar):
            f.write(f"{y1:15.5f} {ybar1:15.5f}\n")
    
    if Lgraph:
        modulename='myplot2D'   ### FOR MLET
        if modulename not in sys.modules:
            import myplot2D
        err, res_max = myplot2D.draw_dots_two(y, y_bar, title, suptitle,Ltwinx=Ltwinx)
    else:
        h_conv = np.array(y_bar) * escale
        y_conv = np.array(y) * escale
        diff =  np.subtract(h_conv,y_conv)
        err = np.sqrt((diff**2).mean())
        res_max = abs(max(diff, key=abs))
        #err = rmse(y, y_bar)*escale
    return err, res_max

def f_write(fname, HL, E_conv, f_conv, err, max_res, job, job_index=None):
    outf = fname.split(".")[0] +'H'+ ''.join(str(x) for x in HL) +'E'+ str(E_conv) + 'F'+ "." + job
    with open(outf, "a") as f:
        if job_index == None:
            f.write("{:5.3f}, {:5.3f}\n".format(err, max_res))
        else:
            f.write("{}: {:5.3f} {:5.3f}\n".format(job_index,err,max_res))
    return 0            
    
def amp_jobs(fdata, job, data_int, amp_pes, HL, E_conv, f_conv, Lgraph, ncore, n_mol, Ltwinx):
    total_images = ase.io.read(fdata, index=':')    # can read extxyz, OUTCAR, 
    images_sets = Images(total_images, nsets=data_int)
    #if not os.path.isfile(amp_pes):

    if re.search("pr", job):
        y=[]
        for mol in total_images:
            y.append(mol.get_potential_energy())
        if fdata.endswith('extxyz'):
            mplot_nvector([],y,fdata.split(".")[0],'sample','E(eV)')
        elif fdata == "OUTCAR":
            mplot_nvector([],y,Xtitle='sample',Ytitle='E(eV)')
    ### JOB == TRAINING
    elif re.search("tr",job):
        if isinstance(data_int, int):
            images = images_sets.get_training_images()
        else:
            d_list =  data_int[:2]
            images = images_sets.get_training_images(d_list=d_list)
        print("data training:total sets %d/%d" % (len(images), len(total_images)))
        calc_train_images(images, HL, E_conv, f_conv, ncore)
        ### test after training:: Do not turn on in qsub
        server =  socket.gethostname()
        if server == 'chi':
            if isinstance(data_int, int):
                images = images_sets.get_test_images()
            else:
                if len(data_int) >= 3:
                    if len(data_int) == 3:
                        d_list = data_int[1:]
                    else:
                        d_list = data_int[2:]
                    images = images_sets.get_test_images(d_list=d_list)
                else:
                    print("There is no test set region in -di ")
            title, suptitle = get_title(fdata, HL, E_conv,f_conv,len(total_images), len(images))
            print("data test:total sets %d/%d" % (len(images), len(total_images)))
            rmserr, max_res = calc_test_images(job, images, amp_pes, title, suptitle,Lgraph,ncore,Ltwinx=Ltwinx)
            f_write(fdata, HL, E_conv, f_conv, rmserr, max_res, job)
    ### JOB == TEST
    elif re.search("te",job):
        if isinstance(data_int, int):
            if data_int == 0:
                images = total_images
            else:
                images = images_sets.get_test_images()
        ### for data interval
        else:
            images = images_sets.get_test_images(d_list=data_int)
            
        title, suptitle = get_title( fdata, HL, E_conv, f_conv, len(total_images), len(images))
        print("data test:total sets %d/%d" % (len(images), len(total_images)))
        rmserr, max_res = calc_test_images(job, images, amp_pes, title, suptitle,Lgraph,ncore,nmol=n_mol,Ltwinx=Ltwinx)
        f_write(fdata, HL, E_conv, f_conv, rmserr, max_res, job)

    elif re.search('md',job):
        # use first geometry 
        atoms = ase.io.read(fdata, index='0')
        run_md(atoms)
    return

def main():
    parser = argparse.ArgumentParser(description='run amp with extxyz, OUTCAR: validation is removed ', prefix_chars='-+/')
    parser.add_argument('-f', '--infile', help='ASE readible file: extxyz, OUTCAR(VASP) ')
    parser.add_argument('-j', '--job', default='tr', choices=['tr','te','pr','pt'], help='job option:"train","test","md","profile","plot"')
    #parser.add_argument('-a', '--all_fig', action="store_true", help='if job==te, include all figures')
    parser.add_argument('-p', '--pot', default="amp.amp", help="input amp potential")
    ### data selection
    data_group = parser.add_argument_group()
    #data_group.add_argument('-d','--data_mine', action="store_true", help='if true, use data interval')
    data_group.add_argument('-di','--data_list', type=int, nargs='+', help='data interval list')
    data_group.add_argument('-dm','--data_mine', type=int, nargs='+', help='data interval list')
    data_group.add_argument('-ds','--dnsets',default=5,type=int, help='num of sets:1 train all sets, otherwise, last set is for test')

    ### others
    parser.add_argument('-nm','--nmol',default=1,type=int,help='num of molecules in the system to normalize error')
    parser.add_argument('-hl', '--hidden_layer', nargs='*', type=int, default=[8,8,8], help='Hidden Layer of lists of integer')
    parser.add_argument('-el', '--e_conv', default=0.001, type=float, help='energy convergence limit')
    parser.add_argument('-fl', '--f_conv', default=0.01, type=float, help='force convergence limit')
    parser.add_argument('-tx', '--twinx', action="store_false", help='turn off to use twinx of two y-axes')
    parser.add_argument('-g', action="store_false", help='if val default is False, otherwise True')
    parser.add_argument('+g', action="store_true", help='if val default is False, otherwise True')
    parser.add_argument('-nc', '--ncore', default=1, type=int, help='number of core needs to be defined')
    args = parser.parse_args()

    if args.data_list:
        data_int = args.data_list
    else:
        data_int = args.dnsets

    pwd = os.getcwd()
    if args.infile:
        fname=args.infile
    else:
        if os.path.isfile('OUTCAR'):
            fname = 'OUTCAR'
        else:
            for f in os.listdir(pwd):
                if f.endswith('extxyz'):
                    fname = f
                    break
        quest = f"file {fname} will be read with job='{args.job}' ?"
        if yes_or_no(quest):
            pass
        else:
            sys.exit(1)

    amp_jobs(fname,args.job,data_int,args.pot,args.hidden_layer,args.e_conv,args.f_conv,args.g,args.ncore,args.nmol,args.twinx)
    return

if __name__ == '__main__':
    main()

