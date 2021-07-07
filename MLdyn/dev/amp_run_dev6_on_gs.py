#!/home/joonho/anaconda3/bin/python

import argparse

import ase.io
from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model.neuralnetwork import NeuralNetwork
from amp.model import LossFunction
from amp.utilities import Annealer

import numpy as np
#from myplot2D import *                 ### FOR MLET 
#import my_chem
#from my_arith import rmse
from my_images import Images
from common import yes_or_no, whereami
from amp_plot import get_title         ### FOR MLET
import re
import sys
import os
import socket

Ldebug = False
amp_amp = [ "amp.amp", "amp-untrained-parameters.amp" ] 


### Train Images 
def calc_train_images(images, HL, Elist, f_conv, f_coeff, ncore, Lload_amp, amp_pot=None):
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
    calc = Amp(descriptor=Gaussian(), model=NeuralNetwork(hiddenlayers=Hidden_Layer), cores=cores)
    ### Global Search in Param Space
    Annealer(calc=calc, images=images, Tmax=20, Tmin=1, steps=4000)

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
    calc.train(images=images, overwrite=True)
    return

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

def calc_test_images(test_images, amp_pes, Ltest_f, title, suptitle,ncore, Lgraph,val_id=None,na_mol=1,Ltwinx=None):
    if amp_pes:
        calc = Amp.load(amp_pes)
    else:
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
        yf=[]
        yf_bar=[]
        f_rmses=[]
    #return         # this runs
    ### as for all the same molecule
    natom = len(test_images[0])
    for mol in test_images:
        y.append(mol.get_potential_energy()/natom*na_mol)
        ### Force
        if Ltest_f:
            y_forces = mol.get_forces()              # 2D list
            if not yf:
                yf.append(y_forces)            # return 2D list of forces
        mol.set_calculator(calc)
        y_bar.append(mol.get_potential_energy()/natom*na_mol)
        if Ltest_f:
            ybar_forces = mol.get_forces()
            if not yf_bar:
                yf_bar.append(ybar_forces)
            yf_arr = np.array(y_forces).flatten()
            yfbar_arr = np.array(ybar_forces).flatten()
            f_rmse = np.sqrt(((yfbar_arr-yf_arr)**2).mean())
            f_rmses.append(f_rmse)
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
            f.write("\tx true\thypo\t\ty true\thypo\t\tz true\thypo\n")
            for f_image, fbar_image in zip(yf, yf_bar):
                for f1, fbar1 in zip(f_image, fbar_image):
                    for i in range(3):
                        f.write(f"{f1[i]:10.5f} {fbar1[i]:10.5f}\t")
                    f.write("\n")
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

def f_write(outf, HL, Elist, f_conv,f_coeff, ntotal, dtype, dlist, err=None, max_res=None, job_index=None):
    with open(outf, "a") as f:
        st = ' '.join(str(x) for x in HL)
        f.write(f"{'Hidden Lay':10}:{st:>10}\n")
        E_conv = Elist[0]
        f.write(f"{'Energy Lim':10}:{E_conv:10g}\n")
        if len(Elist) == 2:
            E_maxres = Elist[1]
            f.write(f"{'Energy max':10}:{E_maxres:10g}\n")
        f.write(f"{'Force  Lim':10}:{f_conv:10g}\n")
        f.write(f"{'Force coef':10}:{f_coeff:10g}\n")
        #print(f"ntotal {ntotal} in {whereami()}")
        f.write(f"{'Total Data':10}:{ntotal:10d}\n")
        f.write(f"{'Data  Type':10}:{dtype:>10}\n")
        if dlist:
            st = ' '.join(str(x) for x in dlist)
            f.write(f"{'Data  list':10}:{st:>10}\n")
        if err:
            if job_index == None:
                f.write("{:5.3f}, {:5.3f}\n".format(err, max_res))
            else:
                f.write(f"{'Error RMSE':10}:{err:10.4f}\n")
                f.write(f"{'MAX Residu':10}:{max_res:10.4f}\n")
    return 0            

def data_selection(total_images, dt, dl, job ):
    
    if dt == 'npart':
        if isinstance(dl, int):
            nset =  dl
        elif isinstance(dl, list):
            nset =  dl[0]
        images_sets     = Images(total_images, dt, nset)
        training_images = images_sets.get_training_images()
        test_images     = images_sets.get_test_images()
        return training_images, test_images    
        
    ### indices for training d[0:2] and test d[len(d)-2:]: not call class Images, job is used 
    elif dt == 'int':
        d_list =  dl[:2]
        training_images = total_images[d_list[0]:d_list[1]]
        if len(dl) >= 3:
            if len(dl) == 3:
                d_list = dl[1:]
            elif len(dl) == 4:
                d_list = dl[2:]
            test_images = total_images[d_list[0]:d_list[1]]
        else:
            test_images=[]
            print("There is no test set region ")
        if job == 'tr':
            return training_images, test_images
        ### one interval will be test region
        elif job == 'te':
            return None, training_images
    ### Division by index: some for training and some for test turn by turn in the file
    elif dt == 'div':
        training_images = []
        test_images = []
        i = 0
        divider     = dl[0]
        tr_remainder = dl[1]
        #    print("Wrong in selection data u. -dt 'div'")
        #    sys.exit(44)
        if len(dl)==3:
            te_remainder = dl[2]
        #te_remain = dl[2]
        for image in total_images:
            if i % divider == tr_remainder:
                training_images.append(image)
                if Ldebug: print(f"{i}-th image in training_images")
            if len(dl) == 3:
                if i % divider == te_remainder:
                    test_images.append(image)
                    if Ldebug: print(f"{i}-th image in test_images")
            i+=1
        if job == 'te':
            test_images=training_images
            training_images=None
        return training_images, test_images
    elif dt == 'pick':
        training_images = []
        test_images = []
        i = 0
        j = 0
        if len(dl) == 2:        # Nontype error for dl, why?
            for image in total_images:
                if i < dl[0]:
                    training_images.append(image)
                    if Ldebug: print(f"{j}-th image in training_images")
                elif i < dl[0]+dl[1]:
                    test_images.append(image)
                    if Ldebug: print(f"{j}-th image in test_images")
                else:
                    training_images.append(image)
                    if Ldebug: print(f"{j}-th image in training_images")
                    i = 0
                i+=1
                j+=1
            return training_images, test_images
        else:
            return None, None
def get_total_image(fdata, ndata):
    ### choose data
    if ndata:
        if len(ndata) == 1:
            st = f':{ndata[0]}'
        elif len(ndata) == 2:
            st = f'{ndata[0]}:{ndata[1]}'
    else:
        st = ':'
    #print(st)
    return ase.io.read(fdata, index=st)    # can read extxyz, OUTCAR, 

def amp_jobs(fdata,job,Ltest_f,amp_inp,Lload_amp,HL,E_conv,f_conv_list,ncore,na_mol,ndata,ntrain,dstype,dslist,Lgraph,Ltwinx):
    ### parameter file
    outf = "hyperparam_" + job + ".dat"     # to write input info
    total_images = get_total_image(fdata,ndata)    # can read extxyz, OUTCAR, 
    print(f"dstype {dstype} dslist {dslist} with data {len(total_images)} in {whereami()}: 1st")
    tr_images, te_images = data_selection(total_images, dstype, dslist, job)
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
    elif re.search("tr",job):
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
        f_write(outf, HL, E_conv, f_conv, f_coeff, ntotal, dstype, dslist)
        calc_train_images(tr_images, HL, E_conv, f_conv, f_coeff, ncore, Lload_amp,  amp_pot=amp_inp)
        ### Test after training:: not run in QSUB
        server =  socket.gethostname()
        if server == 'chi' or server == 'login':
            outf = "runamp_te.dat"     # to write input info
            title, suptitle = get_title(fdata, HL, E_conv,f_conv, f_coeff, ntrain, ntest)
            print("data test:train sets %d/%d" % (ntest, ntrain))
            rmserr, max_res = calc_test_images(te_images,amp_inp,Ltest_f,title, suptitle,ncore,Lgraph,Ltwinx=Ltwinx)
            f_write(outf, HL, E_conv, f_conv, f_coeff, ntotal,dstype, dslist, err=rmserr, max_res=max_res)
    ### JOB == TEST
    elif re.search("te",job):
        title, suptitle = get_title( fdata, HL, E_conv, f_conv,f_coeff, ntrain, ntest)
        print(f"data test {ntest}:{ntrain} train in {whereami()}/job==te")
        rmserr, max_res = calc_test_images(te_images,amp_inp,Ltest_f,title, suptitle,ncore,Lgraph,na_mol=na_mol,Ltwinx=Ltwinx)
        f_write(outf, HL, E_conv, f_conv,f_coeff, ntotal,dstype, dslist, err=rmserr, max_res=max_res)

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
    parser.add_argument('-j', '--job', default='tr', choices=['tr','te','tef','pr','pt','md','chk'], help='job option:"train","test","md","profile","plot","check"')
    parser.add_argument('-tef', '--test_force', action='store_true', help="test force")
    ### Neural Network
    parser.add_argument('-nc', '--ncore', default=1, type=int, help='number of core needs to be defined')
    parser.add_argument('-p', '--pot', help="input amp potential")
    parser.add_argument('-load', '--load_pot', action='store_true', help="load potential")
    parser.add_argument('-nam','--mol_atoms',default=1,type=int,help='num of atoms in molecule in the system to normalize error')
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
        amp_jobs(fname,args.job,args.test_force,args.pot,args.load_pot,args.hidden_layer,args.e_conv,args.f_conv,args.ncore,args.mol_atoms,args.ndata_total,args.ndata_train,args.data_s_type,args.data_s_list,args.g,args.twinx)
    return

if __name__ == '__main__':
    main()

