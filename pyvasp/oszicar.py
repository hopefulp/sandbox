#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os
import sys
import subprocess
import numpy as np
from mplot2D import mplot_nvector as mploty1
from mplot2D import mplot_twinx as mploty2
'''
use free energy (F) for total free energy (E)
    E = F + Ekin + (Skin + Spot)

output
    mdlog.dat   grep T= OSZICAR
    md.dat (default) -of
'''
Ltest = 0
def anal_oszicar(files, outf, ys, iy2s, xlabel, ylabels, title, natom, colors):
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
    tempT = []
    Etot = []
    Efree = []
    E0pot = []
    Ekin = []
    ther_pot = []
    ther_kin = []


    with open("mdlog.dat", 'r') as f:
        for i, line in enumerate(f): 
            #print(line)
            li = line.strip().split()
            t.append(int(li[0]))
            tstep.append(i+1) 
            tempT.append(int(float(li[2])))
            Etot.append(float(li[4]))    ### == E0pot (or Efree) + Ekin + ther_pot + ther_kin
            Efree.append(float(li[6]))
            E0pot.append(float(li[8]))
            Ekin.append(float(li[10]))
            ther_pot.append(float(li[12]))
            ther_kin.append(float(li[14]))
            ###            i-time              temp                 total E             total Efree             E0 pot              E Kin                  Thermo E Kin          Thermo E Pot
            print(f"{int(li[0]):>5} {int(float(li[2])):>5} {float(li[4]):12.5f} {float(li[6]):12.5f} {float(li[8]):12.5f} {float(li[10]):10.5f} {float(li[12]):10.5f} {float(li[14]):10.5f}")
    
    osz = {'T':tempT, 'Etot': Etot, 'Efree':Efree, 'E0pot':E0pot, 'Ekin':Ekin, 'Spot':ther_pot, 'Skin':ther_kin}

    ### E0pot w.r.t. average
    potarr  = np.array(E0pot)
    mean    = np.mean(potarr)
    Epot0   = potarr - mean
    osz['Epot0'] =  Epot0

    ### Ekin per atom
    if natom:
        Ekin_nu = np.array(Ekin)/natom
        osz['Ekin1'] = Ekin_nu
    ### write to file
    fout=open(outf, 'w')
    if Ltest:     diff_E0=[];     diff_Efree=[]
    for i in range(len(t)):
        fout.write(f"{t[i]:>5} {tempT[i]:>5} {Etot[i]:12.5f} {E0pot[i]:12.5f} {Ekin[i]:10.5f}\n")
        ### check energy : Use Ef for Etot ( = EF + Ekin + Spot + Skin)
        #print(f"Etot {Etot[i]:12.5f} E0pot {E0pot[i]:12.5f}  E0pot+Ekin {E0pot[i]+Ekin[i]:12.5f} Eall {E0pot[i]+Ekin[i]+ther_pot[i]+ther_kin[i]:12.5f} Thermo {ther_pot[i]+ther_kin[i]:12.5f}")
        if Ltest:
            diff_E0.append(Etot[i]-(E0pot[i]+Ekin[i]+ther_pot[i]+ther_kin[i]))
            diff_Efree.append(Etot[i]-(Efree[i]+Ekin[i]+ther_pot[i]+ther_kin[i]))
    fout.close()

    ### E0: 0.0179.. Efree: 0.00012
    if Ltest:
        print(f"diff: E0 {np.average(diff_E0):12.5f}, Efree {np.average(diff_Efree):12.5f}")

    ### plot OSZICAR column
    yplots =[]
    ylegends =[]
    all_keys=[]
    if len(ys) == 0:
        print("there should be keys for ys with -y")
        sys.exit(0)
    else:
        #len(ys) == 1:
        #print(f"key: {y}, len {len(y.split())}") # why split?
        #if len(y.split()) == 1:
        ### treat several combinations
        for i, y in enumerate(ys):
            ### treat one combination
            #if '+' in y:
            #keys = y.split('/s+')
            keys = y.split('+')
            print(f" {len(keys)} keys in plot {i+1}")
            all_keys.append(keys)
            ysum = []
            leg=""
            ### combine values in osz.keys()
            for key in keys:
                if key in osz.keys():
                    print(f"add {key} to combination")
                    if leg: leg+="+"
                    leg+=f"{key}"
                    ysum.append(osz[key])
                else:
                    print(f"not {key} in osz keys")
                    sys.exit(11)
            ### sum 2d array
            #yarr1d = np.array(ysum).T.sum(axis=1)
            yarr1d = np.array(ysum).sum(axis=0)
            yplots.append(list(yarr1d))
            ylegends.append(leg)
            test = 0 
            if test == 1:
                ndata = 10
                pivot = len(tstep)-ndata

                for i in range(ndata):
                    print(f"Spot {osz['Spot'][pivot+i]} Skin {osz['Skin'][pivot+i]} Sum {yarr1d[pivot+i]}")
         
    #etotnu = np.array(E0pot) + np.array(Ekin)
    #ther = np.array(ther_pot) + np.array(ther_kin)
    #esum = etotnu + ther
    print(f"dimesions for plot: x {np.array(tstep).shape} y {np.array(yplots).shape}, legend {ylegends}")
    if iy2s:
        mploty2(tstep, yplots, iy2s, legend=ylegends, xlabel=xlabel, ylabel=ylabels, title=title, colors=colors)
    else:
        ylabel = " ".join(ylabels)
        print(f"{ylabel}")
        mploty1(tstep, yplots, legend=ylegends, xlabel=xlabel, ylabel=ylabel, title=title, colors=colors)
    #mplot(tstep, [Etot,etotnu, esum], legend=['Etot','Etotptl','Esum'] )
    ### ys.dim = 1
    #mplot(tstep, T)

    return 0

def main():
    parser = argparse.ArgumentParser(description='VASP AIMD analysis: energy key=T Etot Efree E0pot Ekin Spot Skin')
    parser.add_argument('file', nargs='+', help='read logfile or OSZICAR')
    parser.add_argument('-of', '--outf', default='md.dat', help='save md data')
    #parser.add_argument('-y', '--ys', nargs='+', default=['Etot'], choices=['Etot','Efree','E0pot','Ekin','Spot','Skin'], help='y plots: any combination of Etot Efree E0 Ekin Spot Skin')
    parser.add_argument('-y', '--ys', nargs='+', default=['Etot'], help='y plots: array & combine w. + Etot Efree E0 Ekin Spot Skin')
    parser.add_argument('-xl', '--xlabel', default='t [fs]',  help='label of x-axis')
    parser.add_argument('-yl', '--ylabels', default="E [eV]", nargs='*', help='label of y-axis, list for two y-axes')
    parser.add_argument('-t', '--title', help='title of figure')
    parser.add_argument('-iy2', '--iy2s', nargs='*', type=int, help='y2 column as python index')
    parser.add_argument('-c', '--colors', nargs='*', help='color list')
    parser.add_argument('-na', '--natom', type=int, help='Ekin/natom to get average kinetic energy/Natom')
    args = parser.parse_args()
    
    anal_oszicar(args.file, args.outf, args.ys, args.iy2s, args.xlabel, args.ylabels, args.title, args.natom, args.colors) 

if __name__ == '__main__':
    main()
