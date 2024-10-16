#!/home/joonho/anaconda3/bin/python

'''
written by J. Park 2022/01/13
tdos, pdos for one atom list w.r.t. given Fermi/VBM energy
required mode: lm-decomposed (such as l1m0), more for input atomlist style (include -)

input   atom list which starts from 0
output  pdos atom index which starts from 1
'''
import json
import argparse
import re
import os
import sys
import numpy as np
from mod_doscar import nheadline, obtain_doscar_head, change_Bheadline
from libposcar import obtain_atomlist0
from common     import list2str, whereami
from myplot2D   import mplot_nvector, auto_nvector
from parsing    import convert_2lst2D
from NanoCore import *
def get_max(tdos, pdos):
    if not tdos:
        print(f"max TDOS {max(tdos)}")
    elif not pdos:
        print(f"max PDOS {max(pdos)}")
    else:
        print("Error: no dos")
    return 1
### Nanocore index starts from 1
def extract_doscar_onelist1(doscar, ofile, alist1, idos, Ldoserr, eshift, l, m, bheadline, ngrid):
    '''
    from NanoCore atom starts from 1
    works for one list
    alist1      start from 1
    arr_atom    start from 1
    list        [-1] for TDOS
    '''

    arr_atom = np.array(alist1) # in case atom starts from 0: + 1
    if arr_atom[0] == 0:
        f_pre = 'TDOS'
    else:
        f_pre = 'pdos'

    ### make output file name
    if not ofile:
        if f_pre == 'TDOS':
            ofile = f_pre
            if eshift: f_tdos += eshift[0]
        else:
            alstr = str(arr_atom[0])
            if 1 < len(alist1):
                atom_imax = arr_atom.max()
                atom_imin = arr_atom.min()
                natom_in = arr_atom.size - 2
                if natom_in == 0:
                    alstr = alstr + '_' + str(arr_atom[-1])
                else:
                    if arr_atom.size == atom_imax - atom_imax:
                        alstr = alstr + '-' + str(atom_imax)
                    else:
                        alstr = alstr + '_N' + str(arr_atom.size) + '_' + str(atom_imax)
            ofile = f_pre + 'a' + alstr
            if l: ofile += f'l{l}'
            if m: ofile += f'm{m}'
            if eshift: ofile += eshift[0]
        ofile += '.dat'
    #print(f"{whereami():>15}(): {idos+1}-th output file: {ofile}")
    ### 
    iatom       = -1    # 0 for TDOS and iatom counts from 1 to natom
    iatom_in    =  0
    x_ene       = []    # x needs always to write to pdosa___.dat
    dos         = []    # for TDOS and pdos_sum
    with open(doscar, 'r') as f:
        for i, line in enumerate(f):
            if i < nheadline:
                pass
            ### start new block: headline
            elif line == bheadline:
                iatom += 1                  # 0 for TDOS and iatom counts from 1 to natom
                iene = 0                    # index in energy block
                ### escape after getting all the atoms in arr_atomlist
                if iatom_in == len(arr_atom):
                    break
                ### decide atom box is in or not
                if iatom in arr_atom:
                    Lcheck_atom = True
                    iatom_in += 1
                else:
                    Lcheck_atom = False
                if 0 < iatom and f_pre == 'TDOS': # stop if TDOS
                    break
                #print(f"iatom {iatom} and job {job} in {whereami()}()")
                #if iatom==0:
            ### loop in block: ngrid
            else:
                ### save dos (TDOS or PDOS)
                if Lcheck_atom == True:
                    eles = line.strip().split()
                    ### save x for the first atom in atomlist_array
                    if iatom_in == 1:
                        x_ene.append(float(eles[0]))
                    ### if TDOS  # check spin case
                    if iatom ==0:
                        dos.append(float(eles[1]))
                    ### for pdos:2D # check spin case
                    else:### initialize pdos2d at first
                        if 'pdos2d' not in locals():
                            pdos2d = [ [ 0.0 for x in range(len(eles)-1) ] for y in range(int(ngrid)) ]  # pdos2d[ene][lm]
                            #print(f"{whereami():>15}(): size of 2D pdos {len(pdos2d)} {len(pdos2d[0])}")
                        for i in range(len(eles)-1):
                            pdos2d[iene][i] += float(eles[i+1])
                    iene += 1
                else:
                    continue
    ### energy shift for x_ene: eshift==Evalue|Vvalue
    if eshift:
        ene_shift = float(eshift[1:])
        x_ene = list(np.array(x_ene)-ene_shift)
    if f_pre == 'pdos':
        ### pdos_sum
        dos = list(np.sum(np.array(pdos2d), axis=1))   # projected to atom with sum of all lm's
        #print(f"{whereami():>15}(): size of pdos_sum {np.array(dos).shape}")
    ### remove DOSCAR err at 1st ele
    if Ldoserr:
        del x_ene[0]
        del dos[0]
    
    #print(f"{whereami():>15}(): maxdos {max(dos)} in revised DOSCAR")
    ### writing TDOS[].dat pdosa[].dat
    #fout=open(ofile, 'w')
    #for i in range(len(dos)):
    #    fout.write(f"{x_ene[i]:10.3f}  {dos[i]:10.4f}\n")
    #fout.close()
    #print(f"{whereami():>15}(): write atoms {list(arr_atom)} to {ofile} ")
    
    if f_pre == "TDOS":
        legend = f_pre
    else:
        legend='a'+alstr
    return x_ene, dos, legend

def get_list_z(poscar, dz):
    at = io.read_poscar(poscar)
    atoms = at._atoms
    z_coords=[]; indices=[]

    for atom in atoms:
        if not atom[2] in z_coords: z_coords.append(atom[2])

    zmax = max(z_coords)
    zmin = min(z_coords)
    print(f"z-coord min, max: {zmin} {zmax}")
    zi = zmin - dz/2
    zf = zmax + dz/2

    z_coords.sort()
    for z in z_coords:
        temp = []
        for atom in atoms:
            if abs(z-atom[2]) < 0.01: temp.append(atom.get_serial())
        indices.append(temp)
    #print(f"{indices}")
    print(f"{len(z_coords)} == {len(indices)}")
    return z_coords, indices


def pldos(doscar, poscar, dz, eshift, Lplot, xlabel, ylabel, title):
    ### analyze DOSCAR headline: only one time event
    natom, Emax, Emin, ngrid, Ef, bheadline, Ldoserr, Lspin = obtain_doscar_head(doscar)
    
    ### Use NanoCore
    z_coords, indices = get_list_z(poscar, dz)

    if eshift:
        if re.search('f', eshift, re.I):
            Eshift = 'E' + str(Ef)
            Estr = 'F'
        else:
            Eshift = 'V' + eshift
            Estr = f'V{eshift}'
    else:
        Eshift = None
        Estr = ''
    


    #if Lplot:
    #    ys      = []
    legends = []

    ### loop depending on number of pdos file
    idos = 0
    ### idos to get x, Ldoserr whether remove or not 1st energy
    Z = []
    for inds in indices:
        #x, y          extract_doscar_onelist(doscar,ofile,[-1],idos,Ldoserr, Eshift,l, m, bheadline, ngrid)
        E, dos, legend = extract_doscar_onelist1(doscar,None,inds,idos,Ldoserr, Eshift,0, 0, bheadline, ngrid)
        idos += 1
        Z.append(dos)
        legends.append(legend)

    ### For DOS
    #mplot_nvector(x, ys, xlabel=xlabel, ylabel=ylabel, title=title, legend=legends)
    ### For PLDOS
    absZ = np.abs(Z).T
    print (len(E))
    # convert to log10 values + minimum correction to avoid -INF
    Z = np.log10(absZ + 10**-5)
    E = np.array(E)
    # generate meshgrid
    #X, Y = np.meshgrid(np.array(z_coords), E-Ef)
    X, Y = np.meshgrid(np.array(z_coords), E)

    # customized figure 
    import matplotlib.pyplot as plt
    import pylab as plb
    import matplotlib.ticker as ticker
    fig  = plt.figure(figsize=(10,12))
    ax   = plt.axes()
    #plt.yticks(np.arange(-4,4,1))
    levels = np.linspace(1.01*Z.min(), 0.99*Z.max(), 100)
    cmap=plt.cm.get_cmap("jet")
    cset = plb.contourf(X,Y,Z, levels, cmap=cmap)
    plb.colorbar(cset,ticks=[-4,-3,-2,-1,0,1])
    plt.xticks(np.array(z_coords))
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    #ax.set_ylim(-1.5, 2)
    ax.set_ylim(-0.5, 1.5)
    minor_locator = ticker.AutoMinorLocator(5)
    ax.yaxis.set_minor_locator(minor_locator)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()
    fig.savefig(f'pldos{Estr}.png')


def main():
    parser = argparse.ArgumentParser(description='To obtain PLDOS')
    parser.add_argument('doscar', nargs='?', default='DOSCAR', help='read DOSCAR')
    parser.add_argument('poscar', nargs='?', default="POSCAR", help="input POSCAR")
    parser.add_argument('-dz', '--delta_z', default=0.1, type=float, help='use zmax or delta_z')
    parser.add_argument('-e', '--energy_shift', help='[F|value], F for fermi E shift, V for VBM')
    parser.add_argument('-p', '--plot', action='store_true', help='plot pdos')
    plot = parser.add_argument_group(title='PLOT')
    plot.add_argument('-xl', '--xlabel', default='z (A)', help='xlabel in z-coord')
    plot.add_argument('-yl', '--ylabel', default='E (eV)', help='ylabel for Energy (eV)')
    plot.add_argument('-t', '--title', default='SnO2', help='title for plot')

    args = parser.parse_args()
    
    pldos(args.doscar, args.poscar, args.delta_z, args.energy_shift,args.plot,args.xlabel,args.ylabel,args.title) 

if __name__ == '__main__':
    main()
