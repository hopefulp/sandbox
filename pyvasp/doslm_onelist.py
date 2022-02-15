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
from mod_poscar import obtain_atomlist0
from common     import list2str, whereami
from myplot2D   import mplot_nvector, auto_nvector
from parsing    import convert_2lst_2Dlist

def get_filename(idos, f_pre, arr_atom, eshift, l, m):
    if f_pre == 'TDOS':
        ofile = f_pre
        if eshift: f_tdos += eshift[0]
    else:
        alstr = str(arr_atom[0])
        if 1 < arr_atom.size):
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
        if eshift: ofile += eshift[:]
    ofile += '.dat'
    print(f"{whereami():>15}(): {idos+1}-th output file: {ofile}")
    return ofile, alstr

def extract_doscar_onelist(doscar, ofile, alist0, idos, Ldoserr, eshift, l, m ):
    '''
    works for one list
    alist0      start from 0
    arr_atoms   start from 1
    list        [-1] for TDOS
    '''
    natom, Emax, Emin, ngrid, Ef, bheadline, Ldoserr = obtain_doscar_head(doscar)
    ### change atom index which start from 0
    arr_atoms = np.array(alist0) + 1
    ### tdos or ldos 
    if arr_atoms[0] == 0:
        f_pre = 'TDOS'
    else:
        f_pre = 'ldos'
    ### make filename
    if not ofile:
        ofile = get_filename(idos, f_pre, arr_atoms, eshift, l, m)
    ### 
    iatom       = -1    # 0 for TDOS and iatom counts from 1 to natom
    iatom_in    =  0
    x_ene       = []    # x needs always to write to ldosa___.dat
    dos         = []    # for TDOS and pdos_sum
    pdos2d      = []
    with open(doscar, 'r') as f:
        # nlines = 5 + (1 + nene) * (natom+1) 
        lines = f.readlines()
        
        for iatom in arr_atoms: # iatom is atom index in total system
            istart = nheadline + (1 + nene) * iatom
            for i in range(istart, istart+1+nene):
                eles = lines[i].strip().split()
                ### save x for the first atom in atomlist_array
                if idos == 0:
                    x_ene.append(float(eles[0]))
                ### if TDOS  # update for spin case
                if iatom == 0:
                    dos.append(float(eles[1]))
                ### for pdos:2D # update spin case
                else:### initialize pdos2d at first
                    if not pdos2d:
                        pdos2d = [ [ 0.0 for x in range(len(eles)-1) ] for y in range(int(ngrid)) ]       # pdos2d[ene][lm]
                        print(f"{whereami():>15}(): size of 2D pdos {len(pdos2d)} {len(pdos2d[0])}")
                    for i in range(len(eles)-1):
                        pdos2d[iene][i] += float(eles[i+1])
                iene += 1
    ### energy shift for x_ene: eshift==Evalue|Vvalue
    if eshift:
        ene_shift = float(eshift[1:])
        x_ene = list(np.array(x_ene)-ene_shift)
    if f_pre == 'ldos':
        ### pdos_sum
        dos = list(np.sum(np.array(pdos2d), axis=1))   # projected to atom with sum of all lm's
        print(f"{whereami():>15}(): size of pdos_sum {np.array(dos).shape}")
    ### remove DOSCAR err at 1st ele
    if Ldoserr:
        del x_ene[0]
        del dos[0]
    
    print(f"{whereami():>15}(): maxdos {max(dos)} in revised DOSCAR")
    ### writing TDOS[].dat pdosa[].dat
    fout=open(ofile, 'w')
    for i in range(len(dos)):
        fout.write(f"{x_ene[i]:10.3f}  {dos[i]:10.4f}\n")
    fout.close()
    print(f"{whereami():>15}(): write atoms {list(arr_atom)} to {ofile} ")
    
    if f_pre == "TDOS":
        legend = f_pre
    else:
        legend='a'+alstr
    return x_ene, dos, legend

def extract_doscar(doscar, ofile, alist02d, eshift, l, m, Lplot, xlabel, ylabel, title, colors):
    ### analyze DOSCAR headline: only one time event
    '''
    alist02d    atom list index starts from 0
    '''

    
    if eshift:
        if re.search('f', eshift, re.I):
            Eshift = 'E' + Ef
        else:
            Eshift = 'V' + eshift
    else:
        Eshift = None
    
    print(f"{alist02d}")
    if Lplot:
        ys      = []
        legends = []

### loop depending on number of pdos file
    idos = 0

    for i in range(len(alist02d)):
        x, y, legend = extract_doscar_onelist(doscar,ofile,alist02d[i],idos,Ldoserr,Eshift,l,m,bheadline, ngrid)
        if Lplot:
            ys.append(y)
            legends.append(legend)
    ### plot here
    if Lplot:
        mplot_nvector(x, ys, xlabel=xlabel, ylabel=ylabel, title=title, legend=legends, colors=colors)

def main():
    parser = argparse.ArgumentParser(description='To obtain PLDOS')
    parser.add_argument('doscar', nargs='?', default='DOSCAR', help='read DOSCAR')
    parser.add_argument('poscar', nargs='?', default="POSCAR", help="input POSCAR")
    parser.add_argument('-of', '--ofile', help='output filename')
    ldos = parser.add_mutually_exclusive_group()
    ldos.add_argument('-al','--atom_list0', nargs='*', help="list atoms with num and '-', index from 0 ")
    #ldos.add_argument('-al2', '--atom0_2d', type=json.loads, help='list atoms as 2D list')
    ldos.add_argument('-z', '--zintv', nargs='+', type=float, help='atoms between zmax * zmin')
    parser.add_argument('-ash', '--atom_list0_sh', nargs='*', type=int, help='input shape of atom_list0')
    parser.add_argument('-dz', '--delta_z', default=0.1, type=float, help='use zmax or delta_z')
    parser.add_argument('-loc', '--location', default='in', choices=['in', 'out'], help='outside or inside of zmin')
    parser.add_argument('-as', '--atom_species', nargs='+', default=['O'], help='specify atom species')
    parser.add_argument('-e', '--energy_shift', help='[F|value], F for fermi E shift, V for VBM')
    parser.add_argument('-l', '--ql', type=int, help='angular quantum number')
    parser.add_argument('-m', '--qm', type=int, help='magnetic quantum number')
    parser.add_argument('-p', '--plot', action='store_true', help='plot pdos')
    plot = parser.add_argument_group(title='PLOT')
    plot.add_argument('-xl', '--xlabel', default='E (eV)', help='xlabel for DOS (eV)')
    plot.add_argument('-yl', '--ylabel', default='DOS', help='ylabel for DOS (eV)')
    plot.add_argument('-t', '--title', default='SnO2', help='title for plot')
    plot.add_argument('-c', '--colors', nargs='*', default=['r','g','b','k'], help='colors')

    args = parser.parse_args()
    
    ### obtain atom 2D list
    if args.atom_list0:
        if args.atom_list0_sh:
            alist0 = convert_2lst_2Dlist(args.atom_list0, args.atom_list0_sh)
        else:
            ### convert 1D to 2D list
            alist0=[]
            alist0.append(list(map(int, args.atom_list0)))
    else:
        if not args.zintv:
            print("Error: input z values for zmin and zmax (or dz y d0.01) ")
            sys.exit(1)
        else:
            zaxis=[]
            zaxis.append(args.zintv[0]-args.delta_z)
            zaxis.append(args.zintv[-1]+args.delta_z)
            if args.location == 'in':
                print(f"{whereami():>15}(): use z-axis between {zaxis[0]} and {zaxis[1]}")
            else:
                print(f"{whereami():>15}(): use z-axis below {zaxis[0]} and upper {zaxis[1]}")
            alist0 = obtain_atomlist0(zaxis, args.poscar, args.atom_species, args.location)
    
    extract_doscar(args.doscar, args.ofile, alist0, args.energy_shift,args.ql,args.qm,args.plot,args.xlabel,args.ylabel,args.title,args.colors) 

if __name__ == '__main__':
    main()
