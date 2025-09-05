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
from libposcar import obtain_atomilist0_z, obtain_atomilist0_kind
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



def plot_doscar(doscar, ofile, alist02d, eshift, l, m, Lplot, plot_dict, Lvertical):
    '''
    alist02d    groups of atom list in 2D list: [[3,4,5],[10,11]]
                0   stands for atom index starts from 0
                -1  for TDOS
                atom groups determins by atom.shape of -sh
    l           list with the same size of alist02d: [s, p]
    m

    works for one list
    alist0      start from 0
    arr_atoms   start from 1
    list        [-1] for TDOS
    plot_dict   keys    xlabel, ylabel, xlim, title, colors

    dos_data    np.array of shape (nene, 2) = [lene, dos]
        lene    = dos_data[:,0]     list of energy
        dos     = dos_data[:,1]
    pdos2d      gathers dos
    plot x=lene, y=dos in pdos2d
    '''
    Ef, Lspin, Ldoserr = read_doscar(doscar, option='head')
    
    if Lplot:
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
        print(f"Fermi level: {Ef}")

    ### loop for 2d list
    idos = 0
    if Lspin:
        nlm = 18
    else:
        nlm = 9

    ### As for one atomlist for sum

    pdos2d  = []
    if l:
        if len(l) != len(alist02d):
            print(f"length of l and atom list should be same ")
            sys.exit(1)
    ### make all l[i] = 't'
    else:
        l=[]
        for i in range(len(alist02d)):
            l.append('t')
        
    for i, alist0 in enumerate(alist02d):
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
            fname = get_filename(idos, f_pre, arr_atoms, Eshift, l[i], m)
        else:
            fname = ofile[idos]
        ### dos_data is np.array
        dos_data = read_doscar(filename=doscar, atom_indices=arr_atoms, l=l[i])

        ### dos_data returns always lene, dos
        if not 'lene' in locals():          # as for first dos_data
            lene = dos_data[:,0]
        print(f"lene {lene.shape} in {whereami()}() {__name__}")
        ### remove DOSCAR err at 1st ele
        if Ldoserr:
            if i == 0:                      # as for first dos_data
                lene = np.delete(lene, 0, axis=0)
            dos_data = np.delete(dos_data, 0, axis=0)
        print(f"dimensions lene, dos: {len(list(lene))} {len(list(dos_data))} in module {__name__}")

        ### gather only dos data in 2nd column
        pdos2d.append(dos_data[:,1])
        
        ### writing TDOS[].dat pdosa[].dat in the one atom list; to write several dos at one file f.write out of for loop
        fout=open(fname, 'w')
        for i in range(len(lene)):
            fout.write(f"{lene[i]:10.3f}  {dos_data[i,1]:10.4f}\n")
        fout.close()
        print(f"{whereami():>15}(): write atoms {list(arr_atoms)} to {fname} ")
        '''
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
        '''
    vertical=None
    v_legend = None
    if eshift:
        ene_shift = float(Eshift[1:])
        lene = list(np.array(lene)-ene_shift)
    if Lvertical:
        v_legend = 'Fermi level'
        if eshift: 
            vertical = 0.0
        else:
            vertical = float(Ef)
    print(f"plot vertical line at {vertical}")
    if Lplot:
        ### make legend
        if plot_dict.get('legend'):
            if len(plot_dict['legend']) == len(pdos2d):
                pass
            else:
                print(f"number of legends are different from data: {len(plot_dict['legend'])} {plot_dict['legend']}")
                sys.exit(10)
        elif legends:
            plot_dict['legend'] = legends
        else:
            ### make legends
            for i in range(len(pdos2d)):
                legends.append(f"Model-{i}")
            plot_dict['legend'] = legends
            
    ### plot here
    
    if Lplot:
        ### before plot
        print(f"Before plot: ene dim {len(lene)}, pdos2d {np.array(pdos2d).shape}  in function {whereami()}() in module {__name__}")
        mplot_nvector(lene, pdos2d, plot_dict=plot_dict, Lsave=True, vertical=vertical, v_legend=v_legend)
    return 0

def main():
    parser = argparse.ArgumentParser(description='To obtain PLDOS')
    parser.add_argument('doscar', nargs='?', default='DOSCAR', help='read DOSCAR')
    parser.add_argument('poscar', nargs='?', default="POSCAR", help="input POSCAR")
    parser.add_argument('-of', '--ofile', help='output filename')
    ldos = parser.add_mutually_exclusive_group()
    ldos.add_argument('-al','--atom_list0', nargs='+', default=":", help="atom list from 0, use '-', -1 for Tdos ")
    ldos.add_argument('-ak','--atom_kinds', nargs='+', help='atom kinds if plot all the same atoms'), 
    #ldos.add_argument('-al2', '--atom0_2d', type=json.loads, help='list atoms as 2D list')
    ldos.add_argument('-z', '--zintv', nargs='+', type=float, help='atoms between zmax * zmin')
    parser.add_argument('-ash', '--atom_list0_sh', nargs='*', type=int, help='make atoms groups from atom_list0')
    parser.add_argument('-dz', '--delta_z', default=0.1, type=float, help='use zmax or delta_z')
    parser.add_argument('-loc', '--location', default='in', choices=['in', 'out'], help='outside or inside of zmin')
    parser.add_argument('-as', '--atom_species', nargs='+', default=['O'], help='specify atom species')
    parser.add_argument('-e', '--energy_shift', help='[F|value], -eF[f], -e-3.5 for E such as VBM')
    parser.add_argument('-l', '--ql', nargs='*', help='angular quantum number')
    parser.add_argument('-m', '--qm', type=int, nargs='*', help='magnetic quantum number')
    plot = parser.add_argument_group(title='PLOT')
    plot.add_argument('-p', '--plot', action='store_true', help='plot pdos')
    plot.add_argument('-v', '--vertical', action='store_true', help='add vertical line at Fermi level')
    plot.add_argument('-xl', '--xlabel', default='E (eV)', help='xlabel for DOS (eV)')
    plot.add_argument('-yl', '--ylabel', default='DOS', help='ylabel for DOS (eV)')
    plot.add_argument('-xi', '--xlim', nargs=2, type=float, help='xrange xmin, xmax')
    plot.add_argument('-yi', '--ylim', nargs=2, type=float, help='yrange ymin, ymax')
    plot.add_argument('-t', '--title', default='MoS2-NH$_{x}$', help='title for plot')
    plot.add_argument('-c', '--colors', nargs='*', default=['r','g','b','k'], help='colors')
    plot.add_argument('-lg', '--legend', nargs='*', help='input the same number of legends with plot')
    parser.add_argument('-u', '--usage', action='store_true', help='prints usage and exit')
    args = parser.parse_args()

    if args.usage:
        print("Plot DOSCAR::\
              \n\tdoslm.py -p -eF -xi -3.0 3.0 -al 136-150 159-162 -ash 15 4 -lg I Sn -t 'TEA$_{2}$SnI$_{4}$-V${_I}$’\
              \n\t\
              ")
        sys.exit(0)

    ### obtain atom 2D list
    if args.zintv:
        zaxis=[]
        zaxis.append(args.zintv[0]-args.delta_z)
        zaxis.append(args.zintv[-1]+args.delta_z)
        if args.location == 'in':
            print(f"{whereami():>15}(): use z-axis between {zaxis[0]} and {zaxis[1]}")
        else:
            print(f"{whereami():>15}(): use z-axis below {zaxis[0]} and upper {zaxis[1]}")
        alist0 = obtain_atomilist0_z(zaxis, args.poscar, args.atom_species, args.location)
    else:
        ### if there is shape, there must be list
        if args.atom_list0_sh:
            alist0 = convert_2lst2D(args.atom_list0, args.atom_list0_sh)
        ### Make shape
        elif args.atom_kinds:
            d1, alist0 = obtain_atomilist0_kind(args.poscar, args.atom_kinds)
            print(f"1d list {d1} 2d list {alist0}")
        else:
            alist0=[]
            ### make 1D list if no shape
            if args.atom_list0:
                alist0.append(list(map(int, args.atom_list0)))
    ### gather plot dict
    plot_dict={}
    if args.xlabel: plot_dict['xlabel'] = args.xlabel
    if args.ylabel: plot_dict['ylabel'] = args.ylabel
    if args.xlim:   plot_dict['xlim']   = args.xlim
    if args.ylim:   plot_dict['ylim']   = args.ylim
    if args.title:  plot_dict['title']  = args.title
    if args.colors: plot_dict['colors'] = args.colors
    if args.legend: plot_dict['legend'] = args.legend
    ####                  1            2        3             4           5         6      7         8   
    plot_doscar(args.doscar, args.ofile, alist0, args.energy_shift,args.ql,args.qm,args.plot,plot_dict,args.vertical) 

if __name__ == '__main__':
    main()
