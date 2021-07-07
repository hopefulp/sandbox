#!/home/joonho/anaconda3/bin/python
'''
    2020.08.25 GA(genetic algorithm) was encoded by job='trga', 'tega', if 'ga' in job, turn on Lga
    2020.08.25 Ltest-force is deprecated. if there is force training, make a force test
    2020.11.13 find amp-pot in the directory even though -p None
'''

import argparse
import numpy as np
#import my_chem
from my_images import Images
from common import yes_or_no
from ampplot_test import get_title         ### FOR MLET
import re
import sys
import os
import socket
import amp_util
import amp_descriptor as my_des
from common import whereami
import amp_ini 
from myplot2D import mplot_nvector
from amp_datprocess import outcar_trim

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
def calc_train_images(images, des_obj, HL, Elist, flist, ncore, amp_pot=None, max_iter=10000) :
    Hidden_Layer=tuple(HL)
    if Lprint: print("Hidden Layer: {}".format(Hidden_Layer))
    E_conv, E_maxresid = amp_util.decom_ef(Elist)
    if Lprint: 
        print("Energy convergence: {}".format(E_conv))
        print("Energy maxresidue: {}".format(E_maxresid))
    cores={'localhost':ncore}   # 'localhost' depress SSH, communication between nodes
    ### Detect amp.pot
    ### if amp.pot: load
    calc = amp_util.read_amppot(pot=amp_pot, mvpot=True)
    ### descriptor checking?
    ### if load amp.pot, is this necessary?
    if calc == None:
        if des_obj.name == 'gs':
            ### from amp.descriptor.gaussian2 import Gaussian        # original gaussian, modified gaussian2
            gs = des_obj.make_Gs(images[0]) # images[0] to obtain atom symbols
            #print(gs, f"in {whereami()} of {__name__}")
            #calc = Amp(descriptor=Gaussian(), model=NeuralNetwork(hiddenlayers=Hidden_Layer), cores=cores)
            calc = Amp(descriptor=Gaussian(Gs=gs,cutoff=Cosine(des_obj.cutoff),fortran=True), model=NeuralNetwork(hiddenlayers=Hidden_Layer), cores=cores)
            #calc = Amp(descriptor=my_gauss.Gaussian(Gs=gs,cutoff=Cosine(des_obj.cutoff),fortran=False), model=NeuralNetwork(hiddenlayers=Hidden_Layer), cores=cores)
        elif des_obj.name == 'zn':
            from amp.descriptor.zernike import Zernike
            calc = Amp(descriptor=Zernike(), model=NeuralNetwork(hiddenlayers=Hidden_Layer), cores=cores)
        elif des_obj.name == 'bs':
            from amp.descriptor.bispectrum import Bispectrum
            calc = Amp(descriptor=Bispectrum(), model=NeuralNetwork(hiddenlayers=Hidden_Layer), cores=cores)
        ### Global Search in Param Space
        Annealer(calc=calc, images=images, Tmax=20, Tmin=1, steps=4000)
    ### set convergence for LossFunction
    convergence={}
    if E_conv:
        convergence['energy_rmse'] = E_conv
    if E_maxresid:
        convergence['energy_maxresid'] = E_maxresid
    if flist:
        f_conv, f_coeff = amp_util.decom_ef(flist)
        convergence['force_rmse'] = f_conv      # if it is, force_coefficient turns on
        if not f_coeff:
            f_coeff = 0.04                      # default of amp
    else:
        f_coeff = None
        if Lprint:
            print("Energy only training")
            print(f"force limit: {f_conv:5.2f}, force coeff: {f_coeff:5.2f}")
    calc.model.lossfunction = LossFunction(convergence=convergence, force_coefficient=f_coeff)  # setting
    ### this is always working for max_iteration
    regressor = Regressor(optimizer='BFGS', max_iterations=max_iter)    #'L-BFGS-B'
    calc.model.regressor = regressor
    ### to check train finished
    if os.path.isfile(amp_ini.amptrain_finish):
        os.system(f"rm {amp_ini.amptrain_finish}")
    calc.train(images=images, overwrite=True)
    #calc.train(images=images, max_iterations=max_iter, overwrite=True) # not working in case Annealing runs
    ### Note: when train fails, it might stop running this script here making "amp-untrained-parameters.amp"
    ### Leave message for finishing training
    os.system(f"touch {amp_ini.amptrain_finish}")
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
        modulename='myplot2D'   ### FOR MLET
        if modulename not in sys.modules:
            import myplot2D
        e_rmse, e_maxres = myplot2D.draw_amp_twinx(y, y_bar, title, suptitle, natom=natom, Ltwinx=Ltwinx,Ldiff=True)
    return 0

### Amp job 3: MD
def amp_md(atoms, nstep, dt, amp_pot):
    from ase import units
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.md import VelocityVerlet

    traj = ase.io.Trajectory("traj.traj", 'w')
    calc = amp_util.read_amppot(pot=amp_pot)

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

def amp_jobs(fdata,job,descriptor,amp_inp,HL,E_conv,f_conv_list,max_iter,ncore,na_mol,ndata,ntrain,dstype,dslist,Lgraph,Ltwinx):
    total_images = amp_util.get_total_image(fdata,ndata)    # can read extxyz, OUTCAR, 
    #total_images = outcar_trim(total_images)
    #total_images = ase.io.read(fdata,ndata)
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
    ### parameter file
    if Lprint: print(f"dstype {dstype} dslist {dslist} with data {len(total_images)} in {whereami()}: 1st")
    tr_images, te_images = amp_util.data_partition(total_images, dstype, dslist, job)
    if job == 'te':
        ntotal = 0
        ### ntrain is input
    else:
        ntotal = len(total_images)
        ntrain = len(tr_images)     # ntrain is redefined in case of job==train
    ntest  = len(te_images)
    if Lprint: print(f"ntotal {ntotal} ntrain {ntrain} ntest {ntest} in {whereami()}")
    ### Force input in training
    ### AMP running parts
    Lga = False
    if re.search('ga', job):
        Lga = True
    title, suptitle = get_title(fdata, HL, E_conv,f_conv_list, ntrain, ntest, title=None)
    ### JOB == TRAINING
    if re.search("tr",job):
        pwd = os.getcwd()
        if Lprint: print("data training:total sets %d/%d" % (ntrain, ntotal))
        ### from descriptor - keyword argument
        amp_util.f_write(outf, HL, E_conv, f_conv_list, dstype, dslist, descriptor=descriptor)
        calc_train_images(tr_images, descriptor, HL, E_conv, f_conv_list, ncore, amp_pot=amp_inp,max_iter=max_iter)
            
    ### JOB == TEST
    elif re.search("te",job):
        print(f"data test {ntest} images/ per train {ntrain} in {whereami()}:job==te")
        #os.system("touch test_start")
        calc = amp_util.get_amppot(amp_inp)
        amp_util.f_write(outf, HL, E_conv, f_conv_list, dstype, dslist, descriptor=descriptor, calc=calc, images=te_images)
        calc_test_images(te_images,calc,f_conv_list,title, suptitle,ncore,na_mol,Lgraph,outf=outf,Ltwinx=Ltwinx,Lga=Lga)
    return

def find_inputfile(inf, pwd, job):
    if inf:
        return inf
    elif os.path.isfile('OUTCAR'):
        fname = 'OUTCAR'
    else:
        for f in os.listdir(pwd):
            if f.endswith('extxyz'):
                fname = f
                break
        if "fname" not in locals():
            sys.exit(1)
    #quest = f"file {fname} will be read with job='{job}' ?"
    #if yes_or_no(quest):
    #    return fname
    #else:
    return 0
    

def main():
    parser = argparse.ArgumentParser(description='run amp with extxyz, OUTCAR: validation is removed ', prefix_chars='-+/')
    parser.add_argument('-inf', '--infile', default='OUTCAR', help='ASE readible file: extxyz, OUTCAR(VASP) ')
    ### tef for test w. force
    parser.add_argument('-j', '--job', default='tr', choices=['tr','trga','te','tega','pr','pt','md','chk'], help='job option:"train","test","ga" for addtional genetic algorithm, "md","profile","plot","check"')
    ### Descriptor group
    descriptor_group = parser.add_argument_group(title="Descriptor generator")
    descriptor_group.add_argument('-des', '--descriptor', choices=['gs','zn','bs'], help="test new descriptor")
    descriptor_group.add_argument('-pf', '--p_function', default='log10', choices=['log10','powNN'], help="function for parameter interval")
    descriptor_group.add_argument('-pmm', '--param_minmax', nargs=2, default=[0.05, 5.0], type=float, help="min, max for param interval")
    descriptor_group.add_argument('-pn', '--nparam', type=int, default=5, help="num of parameters for descriptor")
    descriptor_group.add_argument('-pmod', '--param_mod', default='orig', choices=['orig','del','couple','mod'], help="modify param to control number of param")
    descriptor_group.add_argument('-c', '--cutoff', type=float, default=6.5, help="initialize cutoff distance here")
    ### Neural Network
    parser.add_argument('-nc', '--ncore', default=1, type=int, help='number of core needs to be defined')
    parser.add_argument('-p', '--pot', help="input amp potential")
    parser.add_argument('-nam','--natoms_molec',type=int, default=3, help='num of atoms in molecule ')
    parser.add_argument('-hl', '--hidden_layer', nargs='*', type=int, default=[8,8,8], help='Hidden Layer of lists of integer')
    parser.add_argument('-el', '--e_conv', nargs='+', default=[0.001,0.003], type=float, help='energy convergence limit')
    parser.add_argument('-fl', '--f_conv', nargs='*', type=float, help="f-convergence('-' for no f train)[f-coeff]")
    #parser.add_argument('-hf', '--hack_force', action='store_true', help="hacking force in detail")
    #parser.add_argument('-fl', '--f_conv', nargs='*', default=[0.1], type=float, help="f-convergence('-' for no f train)[f-coeff]")
    #parser.add_argument('-fc', '--f_coeff', default=0.1, type=float, help='weight of force training')
    parser.add_argument('-tx', '--twinx', action="store_false", help='turn off to use twinx of two y-axes')
    parser.add_argument('-g', action="store_false", help='if val default is False, otherwise True')
    parser.add_argument('+g', action="store_true", help='if val default is False, otherwise True')
    #parser.add_argument('-a', '--all_fig', action="store_true", help='if job==te, include all figures')
    parser.add_argument('-mi', '--max_iteration', type=int, help='stop at max_iteration for comparison')
    ### DATA group
    data_group = parser.add_argument_group(title='Data')
    data_group.add_argument('-nt', '--ndata_total', type=int, nargs='*', help='cut total data: it requires two value for region')
    data_group.add_argument('-ntr', '--ndata_train', type=int, help='in case of te, define ND_train, ND_test is calculated')
    data_group.add_argument('-dtype','--data_s_type',default='pick',choices=['npart','int','div','pick'], help='data selection type: div-divide by dl[0] and remainder dl[1] for train, dl[2] for test ')
    data_group.add_argument('-dl','--data_s_list', type=int, nargs='+', help='data selection list')
    ### MD group
    md_group = parser.add_argument_group(title='MD')
    md_group.add_argument('-i','--index', default=0, help='select start configuration from input file')
    md_group.add_argument('-ns','--nstep', default=100, type=int, help='number of step with dt')
    md_group.add_argument('-dt','--dt', default=1.0, type=float, help='time interval in fs')

    args = parser.parse_args()

    pwd = os.getcwd()
    fname = find_inputfile(args.infile, pwd, args.job)

    if args.job == 'md':
        index = args.index
        atoms = ase.io.read(fname, index=index)      # index = string
        amp_md(atoms, args.nstep, args.dt, args.pot)
    elif args.job == 'chk':
        pass
    ### if not MD, it's training/test
    else:
        des_obj = my_des.GS_param(args.p_function, pmin=args.param_minmax[0], pmax=args.param_minmax[1], pnum=args.nparam, pmod=args.param_mod, cutoff=args.cutoff)
        amp_jobs(fname, args.job, des_obj, args.pot, args.hidden_layer, args.e_conv, args.f_conv, args.max_iteration, 
            args.ncore, args.natoms_molec, args.ndata_total, args.ndata_train, args.data_s_type, args.data_s_list, args.g, args.twinx)
        #amp_jobs(fname,args.job,des_obj,args.pot,args.hidden_layer,args.e_conv,args.f_conv,args.ncore,args.mol_atoms,args.ndata_total,args.ndata_train,args.data_s_type,args.data_s_list,args.g,args.twinx)
    return

if __name__ == '__main__':
    main()

