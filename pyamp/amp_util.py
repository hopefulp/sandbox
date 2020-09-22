"""
    amp module for amp_run.py

    decom_ef([E_conv_list|force_conve_list)
        decompose arg to [E_rmse|F_rmse], [E_maxres|Force_coefficient]
    get_total_image(fname, ndata)
        return total_images
    data_partition(total_image, data_type, data_list, job)
        return training_images, test_images
    get_inputfile(pwd)
        find 'OUTCAR' or *.extxyz in amp running directory
    f_write(outfilename, HL, Elist, f_conv, f_coeff, ntotal_images, data_type, data_list, descriptor_selection, err, max_res, job_index)
        write arguments to "outfilename"
    write_force_rmse:
    write_energy_diff:

"""
import numpy as np
from amp import Amp
import ase.io
import os
import pickle
from common import whereami, fname_decom
import my_stat
import amp_ini

Ldebug = False
Lprint = 0

def decom_ef(list_ef):
    '''
    force list: returns force_convergence limit, force coefficient
                if force_convergence limit <= 0, it returns None
                if force coefficient returns None, it will be default of amp, 0.04
    energy list: returns energy_convergence limit, energy_maxres
    '''
    ### 1st ele is energy(force)_rmse_limit
    if list_ef == None:
        return None, None
    ef_limit = list_ef[0]
    if ef_limit <= 0:
        ef_limit = None
    ### 2nd ele is energy_maxres or force coefficient: default(amp) = 0.04
    ef_2nd = None
    if len(list_ef) == 2:
         ef_2nd = list_ef[1]
    return ef_limit, ef_2nd

def write_force_rmse(fout, f_rmses):
    with open(fout, 'w') as f:
        f.write(f"{'average':^10}{np.array(f_rmses).mean():10.5f}\n")
        for idx, rmse in enumerate(f_rmses):
            f.write(f"{idx:^9d}{rmse:10.5f} \n")
    return 0            

ffmt = '10.4f'
sfmt  = '^8'

def write_energy(fout, y, y_bar, scale=np.power(10,6)):

    f_pre, f_suf = fname_decom(fout)
    pklfile      = f_pre+ '.pkl'
    p = my_stat.stat_2col(y_bar, y, scale=scale)
    with open(pklfile, 'wb') as f:
        pickle.dump(p, f)
    with open(fout, 'w') as f:
        f.write(f"{'MSE':{sfmt}}{p.mse:{ffmt}}\n{'BIAS_SQ':{sfmt}}{p.bias_sq:{ffmt}}\
                \n{'var(X-Y)':{sfmt}}{p.varx_y:{ffmt}}    {'var(true)':{sfmt}}{p.varx:{ffmt}}    {'var(hypo)':{sfmt}}{p.vary:{ffmt}}    {'Cov(X,Y)':{sfmt}}{p.cov:{ffmt}}    r_corr {p.r_corr:7.4f}\n")
        f.write("\ttrue\t\thypothesis\t\tdiff\n")
        for y1, ybar1 in zip(y, y_bar):
            f.write(f"{y1:15.5f} {ybar1:15.5f} {ybar1-y1:15.5f}\n")
    ### return for score in GA            
    return p.rmse, p.maxres


### sum all files
def stat_forces_filesum(yf3d, ybar_f3d, fpre=None):
    file_force  = "test_force_stat.dat"
    file_maxres = "test_fmaxres"
    accfile     = "test_fstat_acc.txt"
    pklfile     = "test_fstat_acc.pkl"
    f = open(file_force, 'w')

    ### variables for accumulation of images
    true_f = []
    hypo_f = []
    ### loop for images
    i_maxres_file = -1
    tmp_maxres = 0
    for i_image, (yf2d, ybar_f2d) in enumerate(zip(yf3d, ybar_f3d)):
        #print(f"write {fout} in {whereami()}")
        yf_1d       = np.array(yf2d).flatten()
        yfbar_1d    = np.array(ybar_f2d).flatten()
        p           = my_stat.stat_2col(yf_1d, yfbar_1d)
        diff        = np.absolute(yfbar_1d-yf_1d)
        ### select maximum residue
        f_maxres        = np.max(diff)
        f_maxres_ind    = np.argmax(diff)
        if tmp_maxres < f_maxres:
            i_maxres_file  = i_image
            i_maxres_ind   = f_maxres_ind
            tmp_maxres = f_maxres
        ### save for cummulating all the images
        true_f.extend(yf_1d)
        hypo_f.extend(yfbar_1d)

        ### Write force file for each image if fpre != None
        if fpre:
            fout = fpre + '_{:03d}'.format(i_image) + '.dat'
            write_force_image(fout, p, yf2d, ybar_f2d, f_maxres, f_maxres_ind)
        ### write stat to 'test_force_stat.dat'
        f.write(f"{'MSE':{sfmt}}{p.mse:{ffmt}}    {'RMSE':{sfmt}}{p.rmse:{ffmt}}\n{'BIAS_SQ':{sfmt}}{p.bias_sq:{ffmt}}\
            \n{'var(X-Y)':{sfmt}}{p.varx_y:{ffmt}}    {'var(true)':{sfmt}}{p.varx:{ffmt}}    {'var(hypo)':{sfmt}}{p.vary:{ffmt}}    {'Cov(X,Y)':{sfmt}}{p.cov:{ffmt}}    r_corr {p.r_corr:7.4f}\n")
    ### end of loop for files

    ### write maxres file to 'test_fmaxres_img000.txt'
    i_file = i_maxres_file
    yf_1d = np.array(yf3d[i_file]).flatten()
    yfbar_1d = np.array(ybar_f3d[i_file]).flatten()
    p = my_stat.stat_2col(yf_1d, yfbar_1d)
    file_maxres += f"_img{i_maxres_file:03d}.dat"
    write_force_image(file_maxres, p, yf3d[i_file], ybar_f3d[i_file], tmp_maxres, i_maxres_ind)
    f.close()

    ### sum all files
    p = my_stat.stat_2col(true_f, hypo_f, scale=1)
    with open(pklfile, 'wb') as f:
        pickle.dump(p, f)
    with open(accfile, 'w') as f:
        f.write(f"{'MSE':{sfmt}}{p.mse:{ffmt}}\n{'BIAS_SQ':{sfmt}}{p.bias_sq:{ffmt}}\
                \n{'var(X-Y)':{sfmt}}{p.varx_y:{ffmt}}    {'var(true)':{sfmt}}{p.varx:{ffmt}}    {'var(hypo)':{sfmt}}{p.vary:{ffmt}}    {'Cov(X,Y)':{sfmt}}{p.cov:{ffmt}}    r_corr {p.r_corr:7.4f}\n")
    os.system(f"touch {amp_ini.ampout_te_f_chk}")
    return p.rmse, p.maxres

def write_force_image(fout, p, yf2d, ybar_f2d, f_maxres, f_maxres_ind):
    with open(fout, 'w') as f:
        iatom = int(f_maxres_ind/3)
        icoord = f_maxres_ind%3
        f.write(f"{'MSE':{sfmt}}{p.mse:{ffmt}}    {'RMSE':{sfmt}}{p.rmse:{ffmt}}\n{'BIAS_SQ':{sfmt}}{p.bias_sq:{ffmt}}\
            \n{'var(X-Y)':{sfmt}}{p.varx_y:{ffmt}}    {'var(true)':{sfmt}}{p.varx:{ffmt}}    {'var(hypo)':{sfmt}}{p.vary:{ffmt}}    {'Cov(X,Y)':{sfmt}}{p.cov:{ffmt}}    r_corr {p.r_corr:7.4f}\n")
        f.write(f"atom id={iatom} {icoord}-th coord; f_maxres={f_maxres:10.5f}\n")
        f.write(f"{'x true':^10} {'hypo':^10}   {'y true':^10} {'hypo':^10}    {'z true':^10} {'hypo':^10} {'true netforce'} {'hypo_force'} {'diff'}\n")
        #for f_image, fbar_image in zip(yf, yf_bar):   # write all the images?
        atom_f=[]
        for f1, fbar1 in zip(yf2d, ybar_f2d):
            arrf1 = np.array(f1)
            arrfbar1 = np.array(fbar1)
            for k in range(3):
                f.write(f"{f1[k]:10.5f} {fbar1[k]:10.5f}\t")    # write atomic forces in fx, fy, fz
            force_true = (np.sum(arrf1*arrf1))**0.5
            force_hypo = (np.sum(arrfbar1*arrfbar1))**0.5
            diff = force_true - force_hypo
            atom_f.append([force_true, force_hypo])
            f.write(f"{force_true:10.5f} {force_hypo:10.5f} {diff:10.5f}")
            f.write("\n")
    return 0

def f_write(outf, HL, Elist, flist, ntotal, dtype, dlist, descriptor=None, err=None, max_res=None, job_index=None):
    with open(outf, "a") as f:
        if descriptor:
            #f.write(f"{'Descriptor':10}: {descriptor}\n")      # for list, it's simple
            f.write(', '.join("%s: %s" % item for item in vars(descriptor).items()))
            f.write("\n")
        st = ' '.join(str(x) for x in HL)
        f.write(f"{'Hidden Lay':10}:{st:>10}\n")
        E_conv, E_maxres = decom_ef(Elist)
        f.write(f"{'Energy Lim':10}:{E_conv:10g}\n")
        if E_maxres:
            f.write(f"{'Energy max':10}:{E_maxres:10g}\n")
        f_conv, f_coeff = decom_ef(flist)
        if f_coeff:
            f.write(f"{'Force coef':10}:{f_coeff:10g}\n")
        if f_conv:
            f.write(f"{'Force  Lim':10}:{f_conv:10g}\n")
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

def list2str(i_image):
    '''
    limage: (list) range of images
    return (str) ase.io.read argument, '100:120'
    '''
    if len(i_image) == 1:
        sindex = "%i:%i" % (i_image[0], i_image[0]+1)                 # get only one image
        nimage = 1
    else:
        sindex = "%i:%i" % (i_image[0], i_image[1])
        nimage = i_image[1] - i_image[0]
    return sindex

def get_amppot(pot = None):
    if pot:
        calc = Amp.load(pot)
    else:
        try:
            calc = Amp.load("amp.amp")
        except FileNotFoundError:
            try:
                calc = Amp.load("amp-untrained-parameters.amp")
            except FileNotFoundError:
                print("Error: amp-pes.amp file does not exist, input amp-pot file by -p")
                return FileNotFoundError
    return calc

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


def get_inputfile(pwd=None):
    if not pwd:
        pwd = os.getcwd() # this is working also
        print(f"pwd = {pwd} in {whereami()}")
    if os.path.isfile('OUTCAR'):
        fname = 'OUTCAR'
    else:
        for f in os.listdir(pwd):
            if f.endswith('extxyz'):
                fname = f
                break
    return fname

Ntest = 100         # N data for test in case job==training
def data_partition(total_images, dt, dl, job):
    if dt == 'npart':
        if isinstance(dl, int):
            nset =  dl
        elif isinstance(dl, list):
            nset =  dl[0]
        images_sets     = Images(total_images, dt, nset)
        training_images = images_sets.get_training_images()
        test_images     = images_sets.get_test_images()
        return training_images, test_images    
        
    ### INTERVAL: indices for training d[0:2] and test d[len(d)-2:]: not call class Images, job is used 
    elif dt == 'int': 
        d_list=[]
        if not dl:
            if job == 'tr':
                print("There is no test set region ")
            if len(training_images) >= 100:
                training_images = total_images[:-100]
                test_images     = total_images[-100:]
            else:
                test_images=training_images
        elif len(dl) == 1:
            d_list.append(0)
            d_list.append(dl[0])
        elif len(dl) > 4:
            print("data list error")
            sys.exit(10)
        else:
            d_list = dl
        training_images = total_images[d_list[0]:d_list[1]]
        test_images     = total_images[d_list[-2]:d_list[-1]]
        print(f"train data {d_list[0]}:{d_list[1]}, test data {d_list[-2]}:{d_list[-1]}")
        return training_images, test_images
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
