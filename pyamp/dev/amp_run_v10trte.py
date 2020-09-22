#!/home/joonho/anaconda3/bin/python
'''
    2020.08.25 GA(genetic algorithm) was encoded by job='trga', 'tega', if 'ga' in job, turn on Lga
    2020.08.25 Ltest-force is deprecated. if there is force training, make a force test
'''

import argparse
import numpy as np
#from myplot2D import *                 ### FOR MLET 
#import my_chem
#from my_arith import rmse
from my_images import Images
from common import yes_or_no
from amp_plot import get_title         ### FOR MLET
import re
import sys
import os
import socket
import amp_util
import my_descriptor as my_des
from common import whereami
import amp_ini 

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
    ### if amp.pot is loaded
    if amp_pot:
        calc = Amp.load(amp_pot)
    ### descriptor checking?
    if des_obj.name == 'gs':
        #from amp.descriptor.gaussian2 import Gaussian        # original gaussian, modified gaussian2
        gs = des_obj.make_Gs(images[0]) # images[0] to obtain atom symbols
        #print(gs, f"in {whereami()} of {__name__}")
        calc = Amp(descriptor=Gaussian(Gs=gs,cutoff=Cosine(des_obj.cutoff),fortran=False), model=NeuralNetwork(hiddenlayers=Hidden_Layer), cores=cores)
        #calc = Amp(descriptor=my_gauss.Gaussian(Gs=gs,cutoff=Cosine(des_obj.cutoff),fortran=False), model=NeuralNetwork(hiddenlayers=Hidden_Layer), cores=cores)
    elif des_obj.name == 'zn':
        from amp.descriptor.zernike import Zernike
        calc = Amp(descriptor=Zernike(), model=NeuralNetwork(hiddenlayers=Hidden_Layer), cores=cores)
    elif des_obj.name == 'bs':
        from amp.descriptor.bispectrum import Bispectrum
        calc = Amp(descriptor=Bispectrum(), model=NeuralNetwork(hiddenlayers=Hidden_Layer), cores=cores)
    ### Global Search in Param Space
    #Annealer(calc=calc, images=images, Tmax=20, Tmin=1, steps=4000)
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
    calc.train(images=images, overwrite=True)
    #calc.train(images=images, max_iterations=max_iter, overwrite=True) # not working in case Annealing runs 
    return
### AMP job 2: Test
def stat_force_image(y_forces, ybar_forces, f_resmax):
    '''
    to hack force to reduce force_maxresid
    '''
    f_resmax_ind = -1
    yf_1d = np.array(y_forces).flatten()
    yfbar_1d = np.array(ybar_forces).flatten()
    f_rmse = np.sqrt(((yfbar_1d-yf_1d)**2).mean())
    f_resmax_new  = np.max(np.absolute(yfbar_1d-yf_1d))
    print(f"{f_resmax} {f_resmax_new} in {whereami()}")
    if f_resmax < f_resmax_new:
        f_resmax_ind    = np.argmax(np.absolute(yfbar_1d-yf_1d))
        f_resmax = f_resmax_new
    return f_rmse, f_resmax, f_resmax_ind

def calc_test_images(images, amp_pes, Lhack_force, title, suptitle,ncore, Lgraph,Lga=False, val_id=None,na_in_mol=1,Ltwinx=None):
    ### in case other name of amp pot is used
    calc = amp_util.get_amppot(amp_pes)
    escale = 1                  # my_chem.ev2kj
    y=[]            # QM  image energy
    y_bar=[]        # AMP image energy
    yf3d=[]         # QM  force in 3D
    yf3d_bar=[]     # AMP force in 3D
    f_rmses=[]
    f_resmax=0
    ### as for all the same molecule
    ### control display unit: E(ev) per atom(amp), mol, system(all the atoms in image)
    natom = len(images[0]) # Atoms == image, len(Atoms) == number of Atom
    #natom=1
    if na_in_mol == None: na_in_mol = 1
    for i, image in enumerate(images):
        ### get QM-pot
        y.append(image.get_potential_energy()/natom*na_in_mol)
        ### get QM-forces in a image
        if Lhack_force:
            y_forces = image.get_forces()       # 2D list
            yf3d.append(y_forces)               # 3D list 
        ### get AMP-pot
        image.set_calculator(calc)
        y_bar.append(image.get_potential_energy()/natom*na_in_mol)
        ### get AMP-forces
        if Lhack_force:
            ybar_forces = image.get_forces()    # 2D
            yf3d_bar.append(ybar_forces)        # 3D
            ### treat 1 image, obtain iimage of max force_maxresid
            f_rmse, f_resmax, f_resmax_ind = stat_force_image(y_forces, ybar_forces, f_resmax)
            f_rmses.append(f_rmse)
            if f_resmax_ind != -1: 
                f_resmax_iimg = i
                f_resmax_index = f_resmax_ind
    ### write F_rmse (all images): amp_ini.ampout_images_frmse = 'f_rmse.dat'
    if Lhack_force:
        amp_util.write_force_rmse("test_force_rmse.dat", f_rmses)
    ### write energy (all images): 
    amp_util.write_energy("test_energy.txt", y, y_bar)
    ### this runs only in master(login) node
    if Lgraph:
        modulename='myplot2D'   ### FOR MLET
        if modulename not in sys.modules:
            import myplot2D
        err, res_max = myplot2D.draw_amp_twinx(y, y_bar, title, suptitle, natom=natom, Ltwinx=Ltwinx,Ldiff=True)
    else:
        ### write maxforce w. image: amp_ini.ampout_image_maxforce = "forces.dat"
        if Lhack_force:
            amp_util.write_fmaxres_image(f"test_fmaxres_img{f_resmax_iimg}.dat", yf3d[f_resmax_iimg], yf3d_bar[f_resmax_iimg], f_resmax, f_resmax_index)
        ### Running part:: get hl directory from images, amp.amp
        hl = calc.model.parameters.hiddenlayers[images[0][0].symbol]  # retuns hl as tuple
        h_conv = np.array(y_bar)                # * escale
        y_conv = np.array(y)                    # * escale
        diff =  np.subtract(h_conv,y_conv)
        err = np.sqrt((diff**2).mean())
        res_max = abs(max(diff, key=abs))
        #err = rmse(y, y_bar)*escale
        if Lga:
            score = 0.9*res_max+0.1*err
            ### amp_ini.ampout_chromosome_fitness = "ga_fit.txt"
            with open(amp_ini.ampout_chromosome_fitness, 'w') as f:
                f.write(' '.join(str(nodes) for nodes in hl))
                f.write(f" {score}\n")
            os.system(f"cat {amp_ini.ampout_chromosome_fitness} >> ../{amp_ini.ampout_onegeneration_fitness}")  # GA checks this file
        
    return err, res_max


### amp job 3: MD
def amp_md(atoms, nstep, dt, amp_pot):
    from ase import units
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.md import VelocityVerlet

    traj = ase.io.Trajectory("traj.traj", 'w')

    if amp_pot:
        calc = Amp.load(amp_pot)
    else:
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
    dyn = VelocityVerlet(atoms, dt=dt * units.fs)
    f = open("md.ene", "w")
    f.write("{:^5s} {:^10s} {:^10s} {:^10s}\n".format("time","Etot","Epot","Ekin"))
    for step in range(nstep):
        pot = atoms.get_potential_energy()  # 
        kin = atoms.get_kinetic_energy()
        tot = pot + kin
        f.write("{:5d}{:10.5f}{:10.5f}{:10.5f}\n".format(step, tot, pot, kin))
        print("{}: Total Energy={}, POT={}, KIN={}".format(step, tot, pot, kin))
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
    outf = "hyperparam_" + job + ".dat"
    ### parameter file
    total_images = amp_util.get_total_image(fdata,ndata)    # can read extxyz, OUTCAR, 
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
    title, suptitle = get_title(fdata, HL, E_conv,f_conv_list, ntrain, ntest)
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
        if re.search('ga', job):
            Lga = True
        pwd = os.getcwd()
        if Lprint: print("data training:total sets %d/%d" % (ntrain, ntotal))
        amp_util.f_write(outf, HL, E_conv, f_conv_list, ntotal, dstype, dslist, descriptor)
        calc_train_images(tr_images, descriptor, HL, E_conv, f_conv_list, ncore, amp_pot=amp_inp,max_iter=max_iter)
        ### TEST after TRAINING:: deprecated, not run in QSUB, 
        ### if train fails to get "amp.amp", the process stops
        #server =  socket.gethostname()
        #outf = "hyperparam_" + job + "_te.dat"
        #if server == 'chi' or server == 'login':
        #    print("data test/train sets %d/%d" % (ntest, ntrain))
        #    rmserr, max_res = calc_test_images(te_images,amp_inp,Lhack_f,title, suptitle,ncore,Lgraph,na_in_mol=na_mol,Ltwinx=Ltwinx)
        #    amp_util.f_write(outf, HL, E_conv, f_conv_list, ntotal,dstype, dslist, descriptor=descriptor, err=rmserr, max_res=max_res)
        ### Do self test: modify to test test-sets
        #else:
        #    rmserr, max_res = calc_test_images(te_images,amp_inp,Lhack_f,title, suptitle,ncore,Lgraph,Lga=Lga,na_in_mol=na_mol,Ltwinx=False)
        #    amp_util.f_write(outf, HL, E_conv, f_conv_list, ntotal,dstype, dslist, descriptor=descriptor, err=rmserr, max_res=max_res)
            
    ### JOB == TEST
    elif re.search("te",job):
        if f_conv_list: Lhack_f = True
        else          : Lhack_f = False
        print(f"data test {ntest} images/ per train {ntrain} in {whereami()}:job==te")
        #os.system("touch test_start")
        rmserr, max_res = calc_test_images(te_images,amp_inp,Lhack_f,title, suptitle,ncore,Lgraph,Lga=Lga,na_in_mol=na_mol,Ltwinx=Ltwinx)
        amp_util.f_write(outf, HL, E_conv, f_conv_list, ntotal,dstype, dslist, descriptor=descriptor,err=rmserr, max_res=max_res)
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
    #parser.add_argument('-j', '--job', default='tr', choices=['tr','trga','te','tega','tef','pr','pt','md','chk'], help='job option:"train","test","ga" for addtional genetic algorithm, "md","profile","plot","check"')
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
    parser.add_argument('-nam','--mol_atoms',type=int,help='num of atoms in molecule: ')
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
            args.ncore, args.mol_atoms, args.ndata_total, args.ndata_train, args.data_s_type, args.data_s_list, args.g, args.twinx)
        #amp_jobs(fname,args.job,des_obj,args.pot,args.hidden_layer,args.e_conv,args.f_conv,args.ncore,args.mol_atoms,args.ndata_total,args.ndata_train,args.data_s_type,args.data_s_list,args.g,args.twinx)
    return

if __name__ == '__main__':
    main()

