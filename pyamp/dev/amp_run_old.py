#!/home/joonho/anaconda3/bin/python

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
amp_amp = [ "amp.amp", "amp-untrained-parameters.amp" ] 

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

### Amp job 1: Train Images 
def calc_train_images(images, des_obj, HL, Elist, f_conv, f_coeff, ncore, Lload_amp, amp_pot=None):
    Hidden_Layer=tuple(HL)
    print("Hidden Layer: {}".format(Hidden_Layer))
    E_conv = Elist[0]
    if len(Elist) == 2:
        E_maxresid = Elist[1]
    else:
        E_maxresid = E_conv*3
    print("Energy convergence: {}".format(E_conv))
    print("Energy maxresidue: {}".format(E_maxresid))
    cores={'localhost':ncore}   # 'localhost' depress SSH, communication between nodes
    ### load "amp.amp"
    if Lload_amp:
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
    Annealer(calc=calc, images=images, Tmax=20, Tmin=1, steps=4000)
    ### set convergence for LossFunction
    convergence={}
    if E_conv:
        convergence['energy_rmse'] = E_conv
    if E_maxresid:
        convergence['energy_maxresid'] = E_maxresid
    if f_conv > 0.0:
        if f_conv:
            convergence['force_rmse'] = f_conv
        if not f_coeff:
            f_coeff = 0.02
    calc.model.lossfunction = LossFunction(convergence=convergence, force_coefficient=f_coeff)  # setting
    max_iter=500
    ### this is always working
    #regressor = Regressor(optimizer='L-BFGS-B', max_iterations=100)
    regressor = Regressor(optimizer='BFGS', max_iterations=max_iter)
    calc.model.regressor = regressor
    calc.train(images=images, overwrite=True)
    ### this is not working with Annealing
    #calc.train(images=images, max_iterations=max_iter, overwrite=True) 
    return
### AMP job 2: Test
def calc_test_images(test_images, amp_pes, Ltest_f, title, suptitle,ncore, Lgraph,val_id=None,na_in_mol=1,Ltwinx=None):
    ### in case other name of amp pot is used
    #if amp_pes:
    #    calc = Amp.load(amp_pes)
    #else:
    try:
        calc = Amp.load("amp.amp")
    except FileNotFoundError:
        try:
            calc = Amp.load("amp-untrained-parameters.amp")
        except FileNotFoundError:
            print("Error: amp-pes.amp file does not exist, input amp-pot file by -p")
            sys.exit(1)
    escale = 1                  # my_chem.ev2kj
    y=[]
    y_bar=[]
    if Ltest_f:
        yf3d=[]
        yf3d_bar=[]
        f_rmses=[]
        f_maxerr=0
    #return         # this runs
    ### as for all the same molecule
    natom = len(test_images[0])
    for i, image in enumerate(test_images):
        ### get QM-pot
        y.append(image.get_potential_energy()/natom*na_in_mol)
        ### get QM-forces in a image
        if Ltest_f:
            y_forces = image.get_forces()       # 2D list
            yf3d.append(y_forces)               # 3D list 
        ### get AMP-pot
        image.set_calculator(calc)
        y_bar.append(image.get_potential_energy()/natom*na_in_mol)
        ### get AMP-forces
        if Ltest_f:
            ybar_forces = image.get_forces()    # 2D
            yf3d_bar.append(ybar_forces)        # 3D
            ### treat 1 image
            yf_1d = np.array(y_forces).flatten()       
            yfbar_1d = np.array(ybar_forces).flatten()
            f_rmse = np.sqrt(((yfbar_1d-yf_1d)**2).mean())
            f_rmses.append(f_rmse)
            ### max force
            f_res_max = np.max(np.absolute(yfbar_1d-yf_1d))
            if(f_maxerr < f_res_max):
                f_maxerr = f_res_max
                f_max_i = np.argmax(np.absolute(yfbar_1d-yf_1d))
                f_max_image_i = i
    if Ltest_f:
        with open('f_rmse.dat', 'w') as f:
            f.write(f"{'average':^10}{np.array(f_rmses).mean():10.5f}\n")
            for idx, rmse in enumerate(f_rmses):
                f.write(f"{idx:^9d}{rmse:10.5f} \n")

    with open("amp_test.txt",'w') as f:
        f.write("\ttrue\t\thypothesis\n")
        for y1, ybar1 in zip(y, y_bar):
            f.write(f"{y1:15.5f} {ybar1:15.5f}\n")
    if Ltest_f:            
        with open("forces.dat",'w') as f:
            iatom = int(f_max_i/3)
            icoord = f_max_i%3
            f.write(f"image={f_max_image_i}; atom id={iatom} {icoord}-th coord; f_max_res={f_maxerr:10.5f}\n")
            f.write("\tx true\t    hypo\t\t y true\t  hypo\t\t\tz true\t   hypo\t\ttrue_force  hypo_force     diff\n")
            #for f_image, fbar_image in zip(yf, yf_bar):   # write all the images?
            atom_f=[]
            for f1, fbar1 in zip(yf3d[f_max_image_i], yf3d_bar[f_max_image_i]):
                arrf1 = np.array(f1)
                arrfbar1 = np.array(fbar1)
                for i in range(3):
                    f.write(f"{f1[i]:10.5f} {fbar1[i]:10.5f}\t")
                force_true = np.sum(arrf1*arrf1)**0.5
                force_hypo = np.sum(arrfbar1*arrfbar1)**0.5
                diff = force_true - force_hypo
                atom_f.append([force_true, force_hypo])
                f.write(f"{force_true:10.5f} {force_hypo:10.5f} {diff:10.5f}")
                f.write("\n")
            arr_aforce = np.array(atom_f)
            print(f"shape of 2d array: {arr_aforce.shape}")
            ave_f = np.sqrt(np.mean((arr_aforce[:,0]-arr_aforce[:,1])**2))
            f.write(f"force rmse = {ave_f:10.5f}")

    if Lgraph:
        modulename='myplot2D'   ### FOR MLET
        if modulename not in sys.modules:
            import myplot2D
        err, res_max = myplot2D.draw_amp_twinx(y, y_bar, title, suptitle, natom=natom, Ltwinx=Ltwinx,Ldiff=True)
    else:
        h_conv = np.array(y_bar)                # * escale
        y_conv = np.array(y)                    # * escale
        diff =  np.subtract(h_conv,y_conv)
        err = np.sqrt((diff**2).mean())
        res_max = abs(max(diff, key=abs))
        #err = rmse(y, y_bar)*escale
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

def amp_jobs(fdata,job,descriptor,Ltest_f,amp_inp,Lload_amp,HL,E_conv,f_conv_list,ncore,na_mol,ndata,ntrain,dstype,dslist,Lgraph,Ltwinx):
    outf = "hyperparam_" + job + ".dat"
    ### parameter file
    total_images = amp_util.get_total_image(fdata,ndata)    # can read extxyz, OUTCAR, 
    print(f"dstype {dstype} dslist {dslist} with data {len(total_images)} in {whereami()}: 1st")
    tr_images, te_images = amp_util.data_partition(total_images, dstype, dslist, job)
    if job == 'te':
        ntotal = 0
        ### ntrain is input
    else:
        ntotal = len(total_images)
        ntrain = len(tr_images)     # ntrain is redefined in case of job==train
    ntest  = len(te_images)
    print(f"ntotal {ntotal} ntrain {ntrain} ntest {ntest} in {whereami()}: 2nd")
    ### Force input in training
    f_conv = f_conv_list[0]
    f_coeff = 0.1 # default
    if len(f_conv_list) == 2:
        f_coeff = f_conv_list[1]
    ### AMP running parts
    if re.search("pr", job):
        y=[]
        for mol in total_images:
            y.append(mol.get_potential_energy())
        if fdata.endswith('extxyz'):
            mplot_nvector([],y,fdata.split(".")[0],'sample','E(eV)')
        elif fdata == "OUTCAR":
            mplot_nvector([],y,Xtitle='sample',Ytitle='E(eV)')
    ### JOB == TRAINING
    elif re.search("tr",job) or re.search('ga', job):
        pwd = os.getcwd()
        ### prevent running tr in tr-directory
        if amp_inp:
            if Lload_amp:
                pass
            else:
                print(f"{amp_inp} can't be overwritten: Stop")
                sys.exit(11)
        else:
            for amp_file in amp_amp:
                if os.path.isfile(amp_file):
                    print(f"There is {amp_file}: Stop")
                    sys.exit(10)
        print("data training:total sets %d/%d" % (ntrain, ntotal))
        amp_util.f_write(outf, HL, E_conv, f_conv, f_coeff, ntotal, dstype, dslist, descriptor)
        calc_train_images(tr_images, descriptor, HL, E_conv, f_conv, f_coeff, ncore, Lload_amp,  amp_pot=amp_inp)
        ### Test after training:: not run in QSUB
        server =  socket.gethostname()
        outf = "runamp_te.dat"     # to write input info
        if server == 'chi' or server == 'login':
            title, suptitle = get_title(fdata, HL, E_conv,f_conv, f_coeff, ntrain, ntest)
            print("data test:train sets %d/%d" % (ntest, ntrain))
            rmserr, max_res = calc_test_images(te_images,amp_inp,Ltest_f,title, suptitle,ncore,Lgraph,na_in_mol=na_mol,Ltwinx=Ltwinx)
            amp_util.f_write(outf, HL, E_conv, f_conv, f_coeff, ntotal,dstype, dslist, descriptor=descriptor, err=rmserr, max_res=max_res)
        ### Do self test
        else:
            rmserr, max_res = calc_test_images(te_images,amp_inp,Ltest_f,title, suptitle,ncore,Lgraph,na_in_mol=na_mol,Ltwinx=Ltwinx)
            amp_util.f_write(outf, HL, E_conv, f_conv, f_coeff, ntotal,dstype, dslist, descriptor=descriptor, err=rmserr, max_res=max_res)
            
    ### JOB == TEST
    elif re.search("te",job):
        title, suptitle = get_title( fdata, HL, E_conv, f_conv,f_coeff, ntrain, ntest)
        print(f"data test {ntest}:{ntrain} train in {whereami()}/job==te")
        rmserr, max_res = calc_test_images(te_images,amp_inp,Ltest_f,title, suptitle,ncore,Lgraph,na_in_mol=na_mol,Ltwinx=Ltwinx)
        amp_util.f_write(outf, HL, E_conv, f_conv,f_coeff, ntotal,dstype, dslist, descriptor=descriptor,err=rmserr, max_res=max_res)
    return

def find_inputfile(pwd, job):
    if os.path.isfile('OUTCAR'):
        fname = 'OUTCAR'
    else:
        for f in os.listdir(pwd):
            if f.endswith('extxyz'):
                fname = f
                break
    quest = f"file {fname} will be read with job='{job}' ?"
    if yes_or_no(quest):
        return fname
    else:
        return 0
    

def main():
    parser = argparse.ArgumentParser(description='run amp with extxyz, OUTCAR: validation is removed ', prefix_chars='-+/')
    parser.add_argument('-f', '--infile', default='OUTCAR', help='ASE readible file: extxyz, OUTCAR(VASP) ')
    ### tef for test w. force
    parser.add_argument('-j', '--job', default='tr', choices=['tr','te','tef','pr','pt','md','chk','ga'], help='job option:"train","test","md","profile","plot","check"')
    parser.add_argument('-tef', '--test_force', action='store_true', help="test force")
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
    parser.add_argument('-lp', '--load_pot', action='store_true', help="load potential")
    parser.add_argument('-nam','--mol_atoms',default=3,type=int,help='num of atoms in molecule: default 3 for H2O')
    parser.add_argument('-hl', '--hidden_layer', nargs='*', type=int, default=[8,8,8], help='Hidden Layer of lists of integer')
    parser.add_argument('-el', '--e_conv', nargs='+', default=[0.001,0.003], type=float, help='energy convergence limit')
    parser.add_argument('-fl', '--f_conv', nargs='+', default=[0.2, 0.04],   type=float, help='force convergence limit, [force coefficient]')
    #parser.add_argument('-fc', '--f_coeff', default=0.1, type=float, help='weight of force training')
    parser.add_argument('-tx', '--twinx', action="store_false", help='turn off to use twinx of two y-axes')
    parser.add_argument('-g', action="store_false", help='if val default is False, otherwise True')
    parser.add_argument('+g', action="store_true", help='if val default is False, otherwise True')
    #parser.add_argument('-a', '--all_fig', action="store_true", help='if job==te, include all figures')
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
    if args.infile:
        fname=args.infile
    else:
        fname=find_inputfile(pwd, args.job)
    if not fname:
        sys.exit(1)

    if args.job == 'md':
        index = args.index
        atoms = ase.io.read(fname, index=index)      # index = string
        amp_md(atoms, args.nstep, args.dt, args.pot)
    elif args.job == 'chk':
        pass
    ### if not MD, it's training/test
    else:
        des_obj = my_des.GS_param(args.p_function, pmin=args.param_minmax[0], pmax=args.param_minmax[1], pnum=args.nparam, pmod=args.param_mod, cutoff=args.cutoff)
        amp_jobs(fname,args.job,des_obj,args.test_force,args.pot,args.load_pot,args.hidden_layer,args.e_conv,args.f_conv,args.ncore,args.mol_atoms,args.ndata_total,args.ndata_train,args.data_s_type,args.data_s_list,args.g,args.twinx)
    return

if __name__ == '__main__':
    main()

