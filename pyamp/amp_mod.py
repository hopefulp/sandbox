"""
    amp module for amp_run.py
    get_total_image(fname, ndata)
        return total_images
    data_partition(total_image, data_type, data_list, job)
        return training_images, test_images
    f_write(outfilename, HL, Elist, f_conv, f_coeff, ntotal_images, data_type, data_list, descriptor_selection, err, max_res, job_index)
        write arguments to "outfilename"
"""
import ase.io

def f_write(outf, HL, Elist, f_conv,f_coeff, ntotal, dtype, dlist, descriptor, err=None, max_res=None, job_index=None):
    with open(outf, "a") as f:
        f.write(f"{'Descriptor':10}: {descriptor}\n")
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
        if len(dl) == 1:
            d_list.append(0)
            d_list.append(dl[0])
        else:
            d_list =  dl[:2]
        training_images = total_images[d_list[0]:d_list[1]]
        if len(dl) >= 3:
            if len(dl) == 3:
                d_list = dl[1:]
            elif len(dl) == 4:
                d_list = dl[2:]
            test_images = total_images[d_list[0]:d_list[1]]
        else:
            print("There is no test set region ")
            if len(training_images) >= 100:
                test_images=training_images[-100:]
            else:
                test_images=training_images

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
