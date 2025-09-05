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
from libdoscar import *
from common     import list2str, whereami
from mplot2D   import mplot_nvector #, auto_nvector
from parsing    import convert_2lst2D

def get_filename(idos, f_pre, arr_atom, eshift, l, m):
    if f_pre == 'TDOS':
        ofile = f_pre
    else:
        ### atom list string for filename
        ### make atom list string: atom index Amin_max_Nnum
        ofile = f'{f_pre}A{str(arr_atom[0])}'
        if arr_atom.size != 1:
            ofile += f'_{str(arr_atom[-1])}'
        ofile += f'_N{arr_atom.size}'
        if l: ofile += f'l{l}'
        if m: ofile += f'm{m}'
    if eshift:
        if re.search('f', eshift, re.I):
            ofile += 'F0'
        else:
            ofile += eshift[:]
    ofile += '.dat'
    print(f"{whereami():>15}(): {idos+1}-th output file: {ofile}")
    return ofile

def make_vert_line(fname, Ene, maxdos=10):
    f0 = open(fname, 'w')
    dE = 0.001
    EL = float(Ene-dE)
    ER = float(Ene+dE)
    f0.write(f"{EL:10.3f}\t{0.00:10.3f}\n")
    f0.write(f"{Ene:10.3f}\t{maxdos:10.3f}\n")
    f0.write(f"{ER:10.3f}\t{0.00:10.3f}\n")
    f0.close()
    return 0



def extract_doscar(doscar, ofile, alist02d, eshift, l, m, Lplot, plot_dict):
    '''
    alist02d    2D list: groups of atom list, [[3,4,5],[10,11]]
    l           list with the same size of alist02d: [s, p]
    m

    works for one list
    alist0      start from 0
    arr_atoms   start from 1
    list        [-1] for TDOS
    plot_dict   keys    xlabel, ylabel, xlim, title, colors
    '''
    natom, Emax, Emin, nene, Ef, bheadline, Ldoserr, Lspin = obtain_doscar_head(doscar)
    
    if Lplot:
        ys      = []
        legends = []

    if eshift:
        if re.search('f', eshift, re.I):
            Eshift = f"F{Ef:5.3f}"
            ### make F0 vertical line at E=0.00
            make_vert_line("E0.dat", 0.00)      # kw maxdos=10 default
        elif re.search('v', eshift, re.I):
            Eshift = f"V{eshift:5.3f}"
        else:
            Eshift = f"E{eshift:5.3f}"
    else:
        Eshift = None
        ### make F0 vertical line at E=F0
        make_vert_line("EF0.dat", Ef)   # kw maxdos=10 default

    ### loop for 2d list
    idos = 0
    x_ene   = []    # x needs always to write to ldosa___.dat
    dos     = []
    iatom   = 0     # index for total atom count for x_ene=[]
    if Lspin:
        nlm = 18
    else:
        nlm = 9
    with open(doscar, 'r') as f:
        # nlines = 5 + (1 + nene) * (natom+1) 
        lines = f.readlines()

        ### As for one atomlist for sum
        for alist0 in alist02d:
            ### change atom index which starts from 0
            arr_atoms = np.array(alist0) + 1
            ### Make filename: prefix 
            if arr_atoms[0] == 0:
                f_pre = 'TDOS'
            elif l:
                f_pre = 'pdos'
            else:
                f_pre = 'ldos'
            ### make filename
            if not ofile:
                fname = get_filename(idos, f_pre, arr_atoms, Eshift, l, m)
            else:
                fname = ofile[idos]
            pdos2d  = []
            for ind_atom in arr_atoms: # ind_atom is atom index in total system
                istart = nheadline + (1 + nene) * ind_atom
                #print(f"istart {istart} nene {nene} ind_atom {ind_atom} in total lines {len(lines)}")
                iene = 0 # loop index: energy block
                for i in range(istart, istart+1+nene):
                    if lines[i] == bheadline:
                        continue
                    eles = lines[i].strip().split()
                    ### save x for for 1st atom
                    if iatom == 0:
                        x_ene.append(float(eles[0]))
                    ### if TDOS  # update for spin case
                    if ind_atom == 0:
                        dos.append(float(eles[1]))
                    ### for pdos:2D # update spin case
                    else:### initialize pdos2d at first
                        if not pdos2d:
                            pdos2d = [ [ 0.0 for x in range(len(eles)-1) ] for y in range(int(nene)) ]       # pdos2d[ene][lm]
                            print(f"{whereami():>15}(): size of 2D pdos {len(pdos2d)} {len(pdos2d[0])}")
                        for i in range(len(eles)-1):
                            pdos2d[iene][i] += float(eles[i+1])
                    iene += 1
                iatom += 1
            ### one atom list is done, plot dos
            if f_pre == 'ldos':
                ### pdos_sum
                dos = list(np.sum(np.array(pdos2d), axis=1))   # projected to atom with sum of all lm's
                print(f"{whereami():>15}(): size of pdos_sum {np.array(dos).shape}")
            if idos == 0:
                if eshift:
                    ene_shift = float(Eshift[1:])
                    x_ene = list(np.array(x_ene)-ene_shift)
            ### remove DOSCAR err at 1st ele
            if Ldoserr:
                if idos == 0:
                    del x_ene[0]
                del dos[0]
            print(f"dimensions x, dos: {len(x_ene)} {len(dos)}")
            ### writing TDOS[].dat pdosa[].dat in the one atom list
            fout=open(fname, 'w')
            for i in range(len(dos)):
                fout.write(f"{x_ene[i]:10.3f}  {dos[i]:10.4f}\n")
            fout.close()
            print(f"{whereami():>15}(): write atoms {list(arr_atoms)} to {fname} ")
            ### energy shift for x_ene: eshift==Evalue|Vvalue
            idos += 1   # the same as ialist
            print(f"{whereami():>15}(): maxdos {max(dos)} in revised DOSCAR")
            print(f"{idos}th-dos was done")
            ### save var for plot in each dos
            if Lplot:
                ys.append(dos)
                if f_pre == "TDOS":
                    legend = f_pre
                else:
                    legend=re.split('\.',fname)[0]
                legends.append(legend)
    if Lplot and not plot_dict.get('legend'):
        plot_dict['legend'] = legends
    ### plot here
    if Lplot:
        mplot_nvector(x_ene, ys, plot_dict=plot_dict, Lsave=True)
    return 0

def main():
    parser = argparse.ArgumentParser(description='To obtain PLDOS')
    parser.add_argument('doscar', nargs='?', default='DOSCAR', help='read DOSCAR')
    parser.add_argument('poscar', nargs='?', default="POSCAR", help="input POSCAR")
    parser.add_argument('-of', '--ofile', help='output filename')
    ldos = parser.add_mutually_exclusive_group()
    ldos.add_argument('-al','--atom_list0', nargs='+', default=":", help="atom list from 0, use '-', -1 for Tdos ")
    #ldos.add_argument('-al2', '--atom0_2d', type=json.loads, help='list atoms as 2D list')
    ldos.add_argument('-z', '--zintv', nargs='+', type=float, help='atoms between zmax * zmin')
    parser.add_argument('-ash', '--atom_list0_sh', nargs='*', type=int, help='input shape of atom_list0')
    parser.add_argument('-dz', '--delta_z', default=0.1, type=float, help='use zmax or delta_z')
    parser.add_argument('-loc', '--location', default='in', choices=['in', 'out'], help='outside or inside of zmin')
    parser.add_argument('-as', '--atom_species', nargs='+', default=['O'], help='specify atom species')
    parser.add_argument('-e', '--energy_shift', help='[F|value], -eF[f], -e-3.5 for E such as VBM')
    parser.add_argument('-l', '--ql', type=int, nargs='*', help='angular quantum number')
    parser.add_argument('-m', '--qm', type=int, nargs='*', help='magnetic quantum number')
    parser.add_argument('-p', '--plot', action='store_true', help='plot pdos')
    plot = parser.add_argument_group(title='PLOT')
    plot.add_argument('-xl', '--xlabel', default='E (eV)', help='xlabel for DOS (eV)')
    plot.add_argument('-yl', '--ylabel', default='DOS', help='ylabel for DOS (eV)')
    plot.add_argument('-xi', '--xlim', nargs=2, type=float, help='xrange xmin, xmax')
    plot.add_argument('-yi', '--ylim', nargs=2, type=float, help='yrange ymin, ymax')
    plot.add_argument('-t', '--title', default='TEA$_{2}$SnI$_{4}$', help='title for plot')
    plot.add_argument('-c', '--colors', nargs='*', default=['r','g','b','k'], help='colors')
    plot.add_argument('-lg', '--legends', nargs='*', help='input the same number of legends with plot')
    parser.add_argument('-u', '--usage', action='store_true', help='prints usage and exit')
    args = parser.parse_args()

    if args.usage:
        print("Plot DOSCAR::\
              \n\tdoslm.py -p -eF -xi -3.0 3.0 -al 136-150 159-162 -ash 15 4 -lg I Sn -t 'TEA$_{2}$SnI$_{4}$-V${_I}$’\
              \n\t\
              ")
        sys.exit(0)

    ### obtain atom 2D list
    if not args.zintv:
        ### if there is shape, there must be list
        if args.atom_list0_sh:
            alist0 = convert_2lst2D(args.atom_list0, args.atom_list0_sh)
        ### No shape
        else:
            alist0=[]
            ### make 1D list if no shape
            if args.atom_list0:
                alist0.append(list(map(int, args.atom_list0)))
    else:
        zaxis=[]
        zaxis.append(args.zintv[0]-args.delta_z)
        zaxis.append(args.zintv[-1]+args.delta_z)
        if args.location == 'in':
            print(f"{whereami():>15}(): use z-axis between {zaxis[0]} and {zaxis[1]}")
        else:
            print(f"{whereami():>15}(): use z-axis below {zaxis[0]} and upper {zaxis[1]}")
        alist0 = obtain_atomlist0(zaxis, args.poscar, args.atom_species, args.location)
    ### gather plot dict
    plot_dict={}
    if args.xlabel: plot_dict['xlabel'] = args.xlabel
    if args.ylabel: plot_dict['ylabel'] = args.ylabel
    if args.xlim:   plot_dict['xlim']   = args.xlim
    if args.ylim:   plot_dict['ylim']   = args.ylim
    if args.title:  plot_dict['title']  = args.title
    if args.colors: plot_dict['colors'] = args.colors
    if args.legends:plot_dict['legend'] = args.legends
    ####                  1            2        3             4           5         6      7         8   
    extract_doscar(args.doscar, args.ofile, alist0, args.energy_shift,args.ql,args.qm,args.plot,plot_dict) 

if __name__ == '__main__':
    main()
