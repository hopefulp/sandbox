#!/home/joonho/anaconda3/bin/python
'''
    2020.08.25 GA(genetic algorithm) was encoded by job='trga', 'tega', if 'ga' in job, turn on Lga
    2020.08.25 Ltest-force is deprecated. if there is force training, make a force test
    2020.11.13 find amp-pot in the directory even though -p None
'''

import argparse
import numpy as np
import re
import sys
import os
import socket

### in case changing amp package name
#import ampm
#sys.modules['amp'] = ampm
from amp import Amp
from amp.model.neuralnetwork import NeuralNetwork
from amp.model import LossFunction
from amp.regression import Regressor
from amp.utilities import Annealer
from amp.descriptor.cutoffs import Cosine
from amp.descriptor.gaussian import Gaussian
import ase
Lprint = 0

### Amp job 1: Train Images 
def calc_train_images(images, HL, Elist=None, flist=None ) :
    Hidden_Layer=tuple(HL)
    if Lprint: print("Hidden Layer: {}".format(Hidden_Layer))
    #cores={'localhost':ncore}   # 'localhost' depress SSH, communication between nodes
    #print(gs, f"in {whereami()} of {__name__}")
    calc = Amp(descriptor=Gaussian(), model=NeuralNetwork(hiddenlayers=Hidden_Layer))
    ### Global Search in Param Space
    Annealer(calc=calc, images=images, Tmax=20, Tmin=1, steps=4000)
    ### set convergence for LossFunction
    convergence={}
    convergence['energy_rmse'] = 0.001
    #convergence['force_rmse'] = 0.2      # if it is, force_coefficient turns on
    
    #calc.model.lossfunction = LossFunction(convergence=convergence, force_coefficient=0.1)  # setting
    calc.model.lossfunction = LossFunction(convergence=convergence)  # setting
    ### this is always working for max_iteration
    regressor = Regressor(optimizer='BFGS')    #'L-BFGS-B'
    calc.model.regressor = regressor
    calc.train(images=images, overwrite=True)
    ### Note: when train fails, it might stop running this script here making "amp-untrained-parameters.amp"
    ### Leave message for finishing training
    print("Train in finished")
    return 0

### AMP job 2: Test
def calc_test_images(images, calc, f_conv_list, title, suptitle,ncore,na_in_mol,Lgraph, outf='test_result.txt', Ltwinx=None, Lga=False, val_id=None):
    ### in case other name of amp pot is used
    #calc = amp_util.get_amppot(amp_pes)
    Lforce = False
    Lhack_force = False
    if f_conv_list:
        Lforce = True
        Lhack_force = True

    escale = 1                  # my_chem.ev2kj
    y=[]            # QM  image energy
    y_bar=[]        # AMP image energy
    yf3d=[]         # QM  force in 3D
    yf3d_bar=[]     # AMP force in 3D
    ### as for all the same molecule
    ### control display unit: E(ev) per atom(amp), mol, system(all the atoms in image)
    natom = len(images[0]) # Atoms == image, len(Atoms) == number of Atom
    print(f"Lforce = {Lforce}, {Lhack_force}")
    for i, image in enumerate(images):
        ### get QM-pot
        y.append(image.get_potential_energy()/natom*na_in_mol)
        ### get QM-forces in a image
        if Lforce:
            y_forces = image.get_forces()       # 2D list
            yf3d.append(y_forces)               # 3D list 
        ### get AMP-pot
        image.set_calculator(calc)
        y_bar.append(image.get_potential_energy()/natom*na_in_mol)
        ### get AMP-forces
        if Lforce:
            ybar_forces = image.get_forces()    # 2D
            yf3d_bar.append(ybar_forces)        # 3D
    ### write energy (all images): scale np.square(meV)
    e_rmse, e_maxres = amp_util.write_energy(amp_ini.ampout_te_e_chk, y, y_bar, scale=np.power(10,3))
    if Lforce:
        ### write force file of averate and of each image if fpre
        f_rmse, f_maxres = amp_util.stat_forces_filesum(yf3d, yf3d_bar,fpre="test_force_img")
        ###### cal Score
        e_aim   = 0.0005
        f_aim   = 0.1       # 0.1 eV/A
        f2_aim  = 0.3       # 3 sigma
        ### higher score is better 
        score = amp_util.get_score(f_rmse, f_maxres, rmse_target=f_aim, maxres_target=f2_aim)
        amp_util.write_result(outf, score, e_rmse, e_maxres, f_err=f_rmse, f_maxres=f_maxres)
    else:
        score = amp_util.get_score_e(e_rmse, e_maxres)
        amp_util.write_result(outf, score, e_rmse, e_maxres)
     
    ### amp_ini.ampout_score = "ga_fit.txt" >> every test reports "score.dat"
    ### get hl directory from images, amp.amp
    hl = calc.model.parameters.hiddenlayers[images[0][0].symbol]  # hl is tuple, 
    with open('score.dat', 'w') as f:
        f.write(' '.join("%2s" % str(nodes) for nodes in hl))
        f.write(f" {score:10.5f}\n")
    if Lga:
        os.system(f"cat {'score.dat'} >> ../{amp_ini.ampout_onegeneration_fitness}")  # GA checks this file
    
    ### this runs only in master(login) node
    if Lgraph:
        #modulename='myplot2D'   ### FOR MLET
        #if modulename not in sys.modules:
        #    import myplot2D
        #e_rmse, e_maxres = myplot2D.draw_amp_twinx(y, y_bar, title, suptitle, natom=natom, Ltwinx=Ltwinx,Ldiff=True)
        e_rmse, e_maxres = draw_amp_twinx(y, y_bar, title, suptitle, natom=natom, Ltwinx=None,Ldiff=True)
    return 0

### Amp job 3: MD
def amp_md(atoms, nstep, dt, amp_pot):
    from ase import units
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.md import VelocityVerlet

    traj = ase.io.Trajectory("traj.traj", 'w')
    calc = Amp.load(amp_pot)

    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    MaxwellBoltzmannDistribution(atoms, 300 * units.kB)
    traj.write(atoms)
    dyn = VelocityVerlet(atoms, dt=dt * units.fs)
    f = open("md.ene", "w")
    f.write(f"{'time':^5s}{'Etot':^15s}{'Epot':^15s}{'Ekin':^10s}\n")
    print(f"   {'time':^5s}{'Etot':^15s}{'Epot':^15s}{'Ekin':^10s}")
    for step in range(nstep):
        pot = atoms.get_potential_energy()  # 
        kin = atoms.get_kinetic_energy()
        tot = pot + kin
        f.write(f"{step:5d}{tot:15.4f}{pot:15.4f}{kin:10.4f}\n")
        print(f"{step:5d}{tot:15.4f}{pot:15.4f}{kin:10.4f}")
        dyn.run(2)
        traj.write(atoms)                   # write kinetic energy, but pot is not seen in ase
    f.close()        


def decom_force(f_conv_list):
    ### 1st ele is force_rmse_limit
    f_limit = f_conv_list[0]
    if f_limit <= 0:
        f_limit = None
    ### 2nd ele is force coefficient: default(amp) = 0.04
    f_coeff = None
    if len(f_conv_list) == 2:
        f_coeff = f_conv_list[1]
    return f_limit, f_coeff

def amp_jobs(fdata,job,HL,E_conv,f_conv_list,Lgraph):
    total_images = ase.io.read(fdata,':')

    ### draw energy profile of OUTCAR
    if re.search("pr", job):
        y=[]
        for image in total_images:
            y.append(image.get_potential_energy())
            #print(*y, sep='\n')
        if fdata.endswith('extxyz'):
            mplot_nvector([],y,fdata.split(".")[0],'sample','E(eV)')
        elif fdata == "OUTCAR":
            #pass
            mplot_nvector([],y,xlabel='index',ylabel='E(eV)')
            #total_images.write('traj.traj')
        return 0
    outf = "hyperparam_" + job + ".dat"
    ### AMP running parts
    ### JOB == TRAINING
    if re.search("tr",job):
        pwd = os.getcwd()
        ### from descriptor - keyword argument
        calc_train_images(total_images, HL, E_conv, f_conv_list)
            
    ### JOB == TEST
    elif re.search("te",job):
        pass
        #calc = amp_util.get_amppot(amp_inp)
        #amp_util.f_write(outf, HL, E_conv, f_conv_list, dstype, dslist, descriptor=descriptor, calc=calc, images=te_images)
        #calc_test_images(te_images,calc,f_conv_list,title, suptitle,ncore,na_mol,Lgraph,outf=outf,Ltwinx=Ltwinx,Lga=Lga)
    return


def main():
    parser = argparse.ArgumentParser(description='run amp with extxyz, OUTCAR: validation is removed ', prefix_chars='-+/')
    parser.add_argument('-inf', '--infile', default='OUTCAR', help='ASE readible file: extxyz, OUTCAR(VASP) ')
    ### tef for test w. force
    parser.add_argument('-j', '--job', default='tr', choices=['tr','trga','te','tega','pr','pt','md','chk'], help='job option:"train","test","ga" for addtional genetic algorithm, "md","profile","plot","check"')
    ### Neural Network
    parser.add_argument('-hl', '--hidden_layer', nargs='*', type=int, default=[8,8,8], help='Hidden Layer of lists of integer')
    parser.add_argument('-el', '--e_conv', nargs='+', default=[0.001,0.003], type=float, help='energy convergence limit')
    parser.add_argument('-fl', '--f_conv', nargs='*', type=float, help="f-convergence('-' for no f train)[f-coeff]")
    parser.add_argument('-g', action="store_false", help='if val default is False, otherwise True')
    parser.add_argument('+g', action="store_true", help='if val default is False, otherwise True')
    #parser.add_argument('-a', '--all_fig', action="store_true", help='if job==te, include all figures')
    ### MD group
    md_group = parser.add_argument_group(title='MD')
    md_group.add_argument('-p', '--pot', help="input amp potential")
    md_group.add_argument('-i','--index', default=0, help='select start configuration from input file')
    md_group.add_argument('-ns','--nstep', default=100, type=int, help='number of step with dt')
    md_group.add_argument('-dt','--dt', default=1.0, type=float, help='time interval in fs')

    args = parser.parse_args()

    pwd = os.getcwd()

    if args.job == 'md':
        index = args.index
        atoms = ase.io.read(args.infile, index=index)      # index = string
        amp_md(atoms, args.nstep, args.dt, args.pot)
    else:
        amp_jobs(args.infile, args.job, args.hidden_layer, args.e_conv, args.f_conv, args.g)
    return

if __name__ == '__main__':
    main()

