#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys, re, string, getopt, optparse, math, time, pprint 
from os import popen
import bgf
import bgftools
import numpy
import operator
import copy
import nutils as nu
import numpy as np
from numpy import arccos
from numpy.linalg import norm
import itertools
import time
import pickle

option = ""; args = ""; bgf_file = ""; trj_file = ""; out_file = ""; 
usage = """
countWaterCNT.py -b bgf_file -t trj_file -f ff_file(s)
"""
version = "130123"

#-----------------
# calculate the distance of specified atoms in every lammps trajectory file
# this is only used during pei analysis by in kim 
# last update: 2012/06/01
#_________________
def get_line(trj_file):
    """
    Open a lazy file stream. This is considered to increase a file reading speed.
    """
    with open(trj_file) as file:
        for i in file:
            yield i

def sortkey(list):
    return list[0];

def getHBond(bgf_file, trj_file, selection = '', silent=False):

    # init
    timestep = 0; l_timestep = []; line = []; n_header = 0;
    t1 = 0; t2 = 0; # clock
    hbond_dat = dict();
    d_crit = 3.5;
    a_crit = 30.0

    myBGF = bgf.BgfFile(bgf_file)
    myTRJ = open(trj_file)
    myTRJ.seek(0)
    myPickle_file = bgf_file[:-4] + ".hbond.pickle"
    myDAT = open(bgf_file[:-4] + ".hbond.count.dat", "w")
    myDAT.write(str(sys.argv) + "\n")

    # how many steps to go?
    wc_trj_file = popen("grep TIMESTEP " + trj_file + " | wc -l ").read()
    n_timestep = int(wc_trj_file.split()[0]);

    print("The trajectory contains " + str(n_timestep) + " timesteps.")

    # Find header of the trajectory file
    while 1:
        templine = myTRJ.readline()
        line.append(templine.strip('\n').strip('ITEM: '))
        n_header += 1
        if "ITEM: ATOMS" in templine:
            break;

    # INITIAL trajectory information
    timestep = int(line[1])
    natoms = int(line[3])
    boxsize = [line[5].split(' ')[0], line[5].split(' ')[1], line[6].split(' ')[0], line[6].split(' ')[1], line[7].split(' ')[0], line[7].split(' ')[1]]
    boxsize = [float(i) for i in boxsize]
    keywords = line[8].strip('ATOMS ')

    # for every shot in the trajectory file update BGF and manipulate
    dumpatom = get_line(trj_file)
    processed_step = 0;

    t1 = t2 = 0; elapsed_time = 0;

    while 1:
        ### Show progress
        t1 = time.time();
        remaining_time = elapsed_time * (n_timestep - processed_step)
        sys.stdout.write('\r' + "Reading timestep.. " + str(timestep) + " (Remaining time: " + "{0:4.1f}".format(remaining_time) + " seconds, " + "{0:4.1f} minutes".format(remaining_time/60) + ")")
        sys.stdout.flush()

        ### Read
        try:
            chunk = [next(dumpatom) for i in range(natoms+n_header)]
        except StopIteration:
            break;

        timestep = int(chunk[1])
        natoms = int(chunk[3])
        boxsize = [chunk[5].split(' ')[0], chunk[5].split(' ')[1], chunk[6].split(' ')[0], chunk[6].split(' ')[1], chunk[7].split(' ')[0], chunk[7].split(' ')[1]]; 
        boxsize = [float(i) for i in boxsize]; boxsize = [(boxsize[1] - boxsize[0]), (boxsize[3] - boxsize[2]), (boxsize[5] - boxsize[4])]
        keywords = chunk[8].split('ATOMS ')[1].strip('\n').split(' ')


        ### update myBGF with trajectory information ###
        natom_bgf = len(myBGF.a)    # number of atoms in BGF file
    
        if not natom_bgf == natoms:
            nu.die("Number of atoms in trajectory file does not match with BGF file.")
    
        mode = ""
        if 'xs' in keywords or 'ys' in keywords or 'zs' in keywords:
            mode = 'scaled'
        elif 'x' in keywords or 'y' in keywords or 'z' in keywords:
            mode = 'normal'
        elif 'xu' in keywords or 'yu' in keywords or 'zu' in keywords:
            mode = 'unwrapped'
    
        # actual coordinate
        coordinfo = chunk[9:]
    
        # assume that coordinfo is similar to ['id', 'type', 'xs', 'ys', 'zs', 'ix', 'iy', 'iz']
        for atomline in coordinfo:
            atomcoord = atomline.split(' ')
            atom = myBGF.getAtom(int(atomcoord[0]))
    
            if mode == 'scaled':
                atom.x = float(atomcoord[2]) * boxsize[0]
                atom.y = float(atomcoord[3]) * boxsize[1]
                atom.z = float(atomcoord[4]) * boxsize[2]
    
            elif mode == 'unwrapped':
                atom.x = float(atomcoord[2])
                atom.y = float(atomcoord[3])
                atom.z = float(atomcoord[4])
    
            elif mode == 'normal':
                try:
                    ix_index = keywords.index('ix')
                    iy_index = keywords.index('iy')
                    iz_index = keywords.index('iz')
                except ValueError:
                    nu.warn("No image information no the trajectory file. Will be treated as unwrapped.")
                    atom.x = float(atomcoord[2])
                    atom.y = float(atomcoord[3])
                    atom.z = float(atomcoord[4])
                else:
                    atom.x = (int(atomcoord[ix_index]) * boxsize[0]) + float(atomcoord[2]) 
                    atom.y = (int(atomcoord[iy_index]) * boxsize[1]) + float(atomcoord[3]) 
                    atom.z = (int(atomcoord[iz_index]) * boxsize[2]) + float(atomcoord[4]) 
    
            try:
                for i in range(0, 3):
                    myBGF.CRYSTX[i] = boxsize[i]
            except:
                pass;
                #nu.warn("Crystal information error: is this file not periodic?")
                
        pbc = myBGF.CRYSTX[:3]

        # apply periodic condition
        myBGF = bgftools.periodicMoleculeSort(myBGF, myBGF.CRYSTX, silent=True)

        ### myBGF update complete! ###

        # get donor and acceptor
        # donors: O in ACT
        A = []; D = [];
        for atom in myBGF.a:
            if selection:
                if "O" in atom.ffType and eval(selection):
                    A.append(atom)
                if "O" in atom.ffType and eval(selection):
                    D.append(atom)
            else:
                if "O" in atom.ffType:
                    A.append(atom)
                if "O" in atom.ffType:
                    D.append(atom)

        if not len(A) or not len(D):
            nu.die("There are no atoms which can make H_bond (especially O atoms)!")


        ### find hydrogen bonds
        # for all pairs of OC-HW
        sys.stdout.write('Hbonds.. ')
        sys.stdout.flush()

        hbonds = []; donors = []; acceptors = []; hydrogens = []; angles = []; distances = [];

        for d_atom in D:
            d = np.array([d_atom.x, d_atom.y, d_atom.z])
            for a_atom in A:
                a = np.array([a_atom.x, a_atom.y, a_atom.z])
                dist = nu.dist(d, a)
                #dist = bgf.distance(d_atom, a_atom)
                if 0.001 < dist < d_crit:
                    # check H angle
                    for ano in d_atom.CONECT:
                        h_atom = myBGF.getAtom(ano)
                        h = np.array([h_atom.x, h_atom.y, h_atom.z])
                        u = h - d; v = a - d; theta = np.dot(u, v) / norm(u) / norm(v); theta = np.degrees(arccos(theta))
                        #theta = bgf.angle(h_atom, d_atom, a_atom, radians=False)
                        if theta < a_crit:
                            hbonds.append([d_atom.aNo, a_atom.aNo])
                            donors.append(d_atom.aNo)
                            acceptors.append(a_atom.aNo)
                            hydrogens.append(h_atom.aNo)
                            distances.append(dist)
                            angles.append(theta)
                            #n_hbond += 1
                            #print("%s %s" % (d_atom.rNo, a_atom.rNo))
                    
        hbond_dat[timestep] = hbonds
        sys.stdout.write('Done                                                        ')
        sys.stdout.flush()

        print("")
        for i in zip(donors, acceptors, hydrogens, distances, angles):
            print i
        sys.exit(0)

        myDAT.write("%d %d\n" % (timestep, len(hbonds)))

        # write output
        #myBGF.saveBGF(bgf_file[:-4] + "." + str(timestep) + ".bgf")
        #myBGF2.saveBGF(bgf_file[:-4] + "." + str(timestep) + ".bgf")

        t2 = time.time()    # time mark
        elapsed_time = t2 - t1;

    pickle_f = open(myPickle_file, 'w')
    pickle.dump(hbond_dat, pickle_f)
    print('')

    myDAT.close()

    return 1

    ### end of function


if __name__ == "__main__":
    bgf_file = ""; trj_file = ""; ff_file = ""; condition = '';

    options, args = getopt.getopt(sys.argv[1:], 'hb:t:c:', ['help', 'bgf=', 'trj=', 'condition='])

    if len(sys.argv) < 2:
        print(usage)
        sys.exit(0)

    print "Requested options: " + str(options)

    for option, value in options:
        if option in ('-h', '--help'):
            print(usage)
            sys.exit(0)
        elif option in ('-b', '--bgf'):
            bgf_file = value
        elif option in ('-t', '--trj'):
            trj_file = value
        elif option in ('-f', '--forcefield'):
            ff_file = str(value).strip()
        elif option in ('-c', '--selection'):
            condition = str(value)
        elif option == NULL:
            print(usage)
            sys.exit(0)
    
    # more options: resname CNT, resname solvent, save BGF

    # main call
    getHBond(bgf_file, trj_file, selection=condition)
