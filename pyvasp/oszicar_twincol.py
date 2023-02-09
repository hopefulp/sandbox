#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os
import subprocess
from myplot2D import mplot_nvector as mplot
from myplot2D import mplot_twinx as mplot2y
import numpy as np

def anal_oszicar(files, outf, ys, ys2):
    os.system('rm mdlog.dat')
    for fname in files:
        if os.path.isdir(fname):
            f = fname + '/OSZICAR'
        else:
            f = fname
        st = f"grep T= {f} >> mdlog.dat"
        os.system(st)
        #subprocess.check_output(st, shell = True)

    t = []
    tstep = []
    T = []
    Etot = []
    F = []
    Epot = []
    Ekin = []
    ther_pot = []
    ther_kin = []


    with open("mdlog.dat", 'r') as f:
        for i, line in enumerate(f): 
            #print(line)
            li = line.strip().split()
            t.append(int(li[0]))
            tstep.append(i+1) 
            T.append(int(float(li[2])))
            Etot.append(float(li[4]))
            F.append(float(li[6]))
            Epot.append(float(li[8]))
            Ekin.append(float(li[10]))
            ther_pot.append(float(li[12]))
            ther_kin.append(float(li[14]))
            print(f"{int(li[0]):>5} {int(float(li[2])):>5} {float(li[4]):12.5f} {float(li[6]):12.5f} {float(li[8]):12.5f} {float(li[10]):10.5f} {float(li[12]):10.5f} {float(li[14]):10.5f}")
            
    ### write to file
    fout=open(outf, 'w')
    for i in range(len(t)):
        fout.write(f"{t[i]:>5} {T[i]:>5} {Etot[i]:12.5f} {Epot[i]:12.5f} {Ekin[i]:10.5f}\n")
        ### check energy
        print(f"Etot {Etot[i]:12.5f} EF {F[i]:12.5f}  Epot+Ekin {Epot[i]+Ekin[i]:12.5f} Eall {Epot[i]+Ekin[i]+ther_pot[i]+ther_kin[i]:12.5f}")
    fout.close()

    osz = {'T':T, 'Etot': Etot, 'F':F, 'Epot':Epot, 'Ekin':Ekin, 'Spot':ther_pot, 'Skin':ther_kin}
    yplots =[]
    ylegends =[]
    for y in ys:
        print(f"key left: {y}, len {len(y.split())}")
        if len(y.split()) == 1:
            for key in osz.keys():
                if y.casefold() == key.casefold():
                    yplots.append(osz[key])
                    ylegends.append(key)
        ### sum several terms
        else:
            ysum=[]
            leg=""
            inp_keys = y.split()
            for inp_key in inp_keys:
                for key in osz.keys():
                    if inp_key.casefold() == key.casefold():
                        print(f"add {key} to ys")
                        ysum.append(osz[key])
                        leg+=f"{key}+"
            ### sum 2d array
            yarr1d = np.array(ysum).T.sum(axis=1)
            yplots.append(yarr1d)
            ylegends.append(leg)
    if ys2:
        yplots2 =[]
        ylegends2 =[]
        for y in ys2:
            print(f"key right: {y}, len {len(y.split())}")
            if len(y.split()) == 1:
                for key in osz.keys():
                    if y.casefold() == key.casefold():
                        yplots2.append(osz[key])
                        ylegends2.append(key)
            ### sum several terms
            else:
                ysum=[]
                leg=""
                inp_keys = y.split()
                for inp_key in inp_keys:
                    for key in osz.keys():
                        if inp_key.casefold() == key.casefold():
                            print(f"add {key} to ys")
                            ysum.append(osz[key])
                            leg+=f"{key}+"
                ### sum 2d array
                yarr1d = np.array(ysum).T.sum(axis=1)
                yplots2.append(yarr1d)
                ylegends2.append(leg)
        
            
    #etotnu = np.array(Epot) + np.array(Ekin)
    #ther = np.array(ther_pot) + np.array(ther_kin)
    #esum = etotnu + ther
    if not ys2:
        mplot(tstep, yplots, legend=ylegends)
    else:
        leg = [ylegends, ylegends2]
        mplot2y(tstep, yplots, yplots2, legend=leg)
    #mplot(tstep, [Etot,etotnu, esum], legend=['Etot','Etotptl','Esum'] )
    ### ys.dim = 1
    #mplot(tstep, T)

    return 0

def main():
    parser = argparse.ArgumentParser(description='read vasp logfile')
    parser.add_argument('file', nargs='+', help='read logfile or OSZICAR')
    parser.add_argument('-of', '--outf', default='md.dat', help='save md data')
    parser.add_argument('-y', '--ys', nargs='+', default='F', help='y plots')
    parser.add_argument('-y2', '--ys2', nargs='*', help='y plots')
    args = parser.parse_args()
    
    anal_oszicar(args.file, args.outf, args.ys, args.ys2) 

if __name__ == '__main__':
    main()
