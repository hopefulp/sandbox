#!/home/joonho/anaconda3/bin/python

'''
written by J. Park 2022/01/13
tdos, pdos for one atom list w.r.t. given Fermi/VBM energy
required mode: lm-decomposed (such as l1m0), multiple list input using json

input   atom list which starts from 0
output  pdos atom index which starts from 1
'''



import argparse
import re
import os
import sys
import numpy as np
from mod_doscar import nheadline, obtain_doscar_head, change_Bheadline
from mod_poscar import obtain_atomlist0
from common import list2str, whereami
from myplot2D import mplot_nvector, auto_nvector

dos_err = 100

def get_max(tdos, pdos):
    if not tdos:
        print(f"max TDOS {max(tdos)}")
    elif not pdos:
        print(f"max PDOS {max(pdos)}")
    else:
        print("Error: no dos")
    return 1

def extract_doscar(doscar, ofile, job, alist0, eshift, l, m, Lplot):
    '''
    alist0      start from 0
    arr_atom    start from 1
    '''
    arr_atom = np.array(alist0) + 1
    f_pre = 'pdos'
    f_tdos = 'TDOS'
    f_suff = ''
    if eshift:
        if re.search('f', eshift[0], re.I):
            f_suff='F'
        elif re.search('v', eshift[0], re.I):
            f_suff='V'

    if not ofile:
        alstr = str(arr_atom[0])
        if 1 < len(alist0):
            alstr = alstr + '_' + str(arr_atom[-1])
        ofile = f_pre + 'a' + alstr
        if l: ofile += f'l{l}'
        if m: ofile += f'm{m}'
        if f_suff: ofile += f_suff
        ofile += '.dat'
    if f_suff: f_tdos += f"{f_suff}.dat"
    print(f"{whereami():>15}(): output file: {ofile}")
    ### analyze DOSCAR headline
    natom, Emax, Emin, ngrid, Ef, bheadline, E2nd = obtain_doscar_head(doscar)
    ### Doesnot need to change bheadline in analysis DOSCAR, just skip the 1st energy in abnormality
    totline = nheadline + (natom+1) * (int(ngrid)+1) 
    iatom=-1     # 0 for TDOS and iatom counts from 1 to natom
    x_ene=[]
    with open(doscar, 'r') as f:
        for i, line in enumerate(f):
            if i < nheadline:
                pass
            ### start new block: headline
            elif line == bheadline:
                iatom += 1                  # 0 for TDOS and iatom counts from 1 to natom
                iene = 0                    # index in energy block
                #print(f"iatom {iatom} and job {job} in {whereami()}()")
                if iatom==0:
                    tdos = []
                    tmax = 0.
            ### loop in block: ngrid
            else:
                ### save TDOS or PDOS
                if iatom==0:
                    eles = line.strip().split()
                    x_ene.append(float(eles[0]))
                    tdos.append(float(eles[1]))
                ### skip and write twice each first line in the energy block
                elif iatom in arr_atom:                # arr_atom: 1 ~
                    eles = line.strip().split()
                    ### at first atom, prepare columns depending on spin
                    if 'pdos2d' not in locals():
                        pdos2d = [ [ 0.0 for x in range(len(eles)-1) ] for y in range(int(ngrid)) ]       # pdos2d[ene][lm]
                        print(f"{whereami():>15}(): size of 2D pdos {len(pdos2d)} {len(pdos2d[0])}")
                    for i in range(len(eles)-1):
                        pdos2d[iene][i] += float(eles[i+1])
                    iene += 1
                else:
                    continue
    ### modified energy in case Efermi or EVBM
    if f_suff:
        if f_suff == 'F':
            ene_shift = Ef
        elif f_suff == 'V':
            ene_shift = float(eshift[1])
        x_ene = list(np.array(x_ene)-ene_shift)

    pdos_sum = list(np.sum(np.array(pdos2d), axis=1))   # projected to atom with sum of all lm's
    print(f"{whereami():>15}(): size of pdos_sum {np.array(pdos_sum).shape}")
    ### test max dos error and remove it from all doses: dos is any of tdos|pdos
    if 'tdos' in locals():
        dos = tdos
    else:
        dos = pdos_sum
    maxdos = max(dos)
    imaxdos = dos.index(maxdos)
    print(f"{whereami():>15}(): maxdos {maxdos} in {imaxdos} and at {imaxdos+1} dos {dos[imaxdos+1]}")
    Lrm_max = 0
    if dos[imaxdos+1] * dos_err < maxdos:
        Lrm_max = 1

    ### reset dos by removing max elements
    ### all the variables to be reported arranged here
    if Lrm_max:
        if 'tdos' in locals():
            del tdos[0]
        del x_ene[0]
        del pdos_sum[0]
        ngrid -= 1
    
    print(f"{whereami():>15}(): maxdos {max(dos)} in revised DOSCAR")
    ### writing TDOS.dat
    if 't' in job:
        fout=open(f_tdos, 'w')
        for i in range(int(ngrid)):
            fout.write(f"{x_ene[i]:10.3f}  {tdos[i]:10.2f}\n")
        fout.close()
        print(f"write {f_tdos}")
    ### writing pdosa#_#.dat
    fout=open(ofile, 'w')
    for i in range(int(ngrid)):
        fout.write(f"{x_ene[i]:10.3f}  {pdos_sum[i]:10.4f}\n")
    print(f"write {arr_atom}: {ofile}")
    fout.close()
    ### Plot
    if Lplot:
        ys=[]
        legend=[]
        if 'tdos' in locals():
            ### tdos should be numbers for plot
            ys.append(tdos)
            legend.append('TDOS')
        ys.append(pdos_sum)
        legend.append('a'+alstr)
        ### x should be numbers not str
        mplot_nvector(x_ene, ys, xlabel='E (eV)', ylabel='DOS', legend=legend)
    return 0

def main():
    parser = argparse.ArgumentParser(description='To obtain PLDOS')
    parser.add_argument('doscar', nargs='?', default='DOSCAR', help='read DOSCAR')
    parser.add_argument('poscar', nargs='?', default="POSCAR", help="input POSCAR")
    parser.add_argument('-of', '--ofile', help='output filename')
    parser.add_argument('-j', '--job', default='t', help='l,p,t, LDOS, PDOS, LPDOS, TDOS')
    ldos = parser.add_mutually_exclusive_group()
    ldos.add_argument('-al','--atom_list0', nargs='*', type=int, help="list atoms: index from 0 ")
    ldos.add_argument('-z', '--zintv', nargs='+', type=float, help='atoms between zmax * zmin')
    parser.add_argument('-dz', '--delta_z', default=0.1, type=float, help='use zmax or delta_z')
    parser.add_argument('-loc', '--location', default='in', choices=['in', 'out'], help='outside or inside of zmin')
    parser.add_argument('-as', '--atom_species', nargs='+', default=['O'], help='specify atom species')
    parser.add_argument('-e', '--energy_shift', nargs='*', help='angular quantum number')
    parser.add_argument('-l', '--ql', type=int, help='angular quantum number')
    parser.add_argument('-m', '--qm', type=int, help='magnetic quantum number')
    parser.add_argument('-p', '--plot', action='store_true', help='plot pdos')
    #parser.add_argument('-o', '--old_suffix', default='_o', help='save original DOSCAR to DOSCAR_o')
    args = parser.parse_args()
    
    ### obtain atom list
    if args.atom_list0:
        alist0 = args.atom_list0
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
    
    extract_doscar(args.doscar, args.ofile, args.job, alist0, args.energy_shift, args.ql, args.qm, args.plot) 

if __name__ == '__main__':
    main()
