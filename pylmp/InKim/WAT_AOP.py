#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys, re, string, getopt, optparse, math, time, pprint 
from os import popen
import bgf
import bgftools
import nutils as nu
import numpy as np
import numpy.linalg as la
import scipy.spatial
import periodic_kdtree as pkdtree
import itertools
import time
#import lammpstrj as lt
from lammps.trj import *

option = ""; args = ""; bgf_file = ""; trj_file = ""; out_file = ""; 
usage = """
countWaterCNT.py -b bgf_file -t trj_file -f ff_file(s)
"""
version = "160516"

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


def countWaterCNT(bgf_file, trj_file, d_crit = 3.5, selection='', nsample = 0, silent=False):

    ### init
    # constants
    deg_109p47 = math.radians(109.47)

    # variables
    timestep = 0; l_timestep = []; line = []; n_header = 0;
    avg_angles = []; avg_diheds = []
    bin = np.arange(0.0, 181.0, 1.0)
    t1 = 0; t2 = 0; # clock

    myBGF = bgf.BgfFile(bgf_file)
    myTRJ = open(trj_file)
    myTRJ.seek(0)
    myDAT = open(bgf_file[:-4] + ".AOP.dat", "w")
    myDAT.write("#timestep\tAOP\tF4\n")

    # how many steps to go?
    #mytrj = lt.lammpstrj(trj_file)
    mytrj = Trj(trj_file)
    l_timestep = mytrj.load()
    n_timestep = len(l_timestep)
    l_timestep.sort()
    requested_timesteps = l_timestep[-nsample:]

    print("Using neighboring water distance criteria %4.1f A (should be consistent with the coordination number)" % d_crit)
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
        remaining_time = elapsed_time * (len(requested_timesteps) - processed_step)
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

        if not timestep in requested_timesteps:
            continue

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
                
        # apply periodic condition
        myBGF = bgftools.periodicMoleculeSort(myBGF, myBGF.CRYSTX, silent=True)

        ### myBGF update complete! ###

        # get water OW
        O = []; anos_O = []
        for atom in myBGF.a:
            if selection:
                if "O" in atom.ffType and eval(selection):
                    O.append(atom)
                    anos_O.append(atom.aNo)
            else:
                if "O" in atom.ffType:
                    O.append(atom)
                    anos_O.append(atom.aNo)
                
        if not len(O):
            nu.die("There are no O atoms which satisfies %s!" % selection)


        ### calculate AOP and F4
        sys.stdout.write('AOP..... ')
        sys.stdout.flush()

        coords = []; anos = [];
        for atom in O:
            coords.append([atom.x, atom.y, atom.z])
            anos.append(atom.aNo)

        pbc = np.array(myBGF.CRYSTX[:3])
        coords = np.array(coords)
        tree = pkdtree.PeriodicKDTree(pbc, coords)  # KDtree for distance calculation

        # analyze local water structures for a timestep
        n_neighbors = []; angles = []; diheds = [];
        for atom in O:
            # local variables

            # find neighbors
            neighbors = []; # list of O atoms near centerO
            centerO = np.array([atom.x, atom.y, atom.z])
            d, ndx = tree.query(centerO, k=6)
            index = np.where((d >= 1e-4) & (d <= d_crit))  # discard self
            for i in ndx[index]:
                neighbors.append(myBGF.getAtom(anos[i]))

            # calculate O1-O-O2 angle
            for i, j in itertools.combinations(neighbors, 2):
                if not "O" in i.ffType or not "O" in j.ffType:
                    nu.die("Wrong atom found.")
                x = [i.x, i.y, i.z]; y = [j.x, j.y, j.z]
                pbc_dist = nu.pbc_dist(x, y, pbc)
                dist = nu.dist(x, y)
                if abs(pbc_dist - dist) > 0.1: continue;    # discard atoms over pbc boundaries
                
                # angle
                angle = nu.angle(x, centerO, y, radians=False)
                angles.append(angle)

            n_neighbors.append(len(neighbors))  # number of neighboring water molecules

            # calculate dihedral
            for i in neighbors:
                # find aNos of H located farthest among all combinations
                hO_anos = atom.CONECT; hi_anos = i.CONECT
                max_dist_pair = []; max_dist = 0.0;
                for k, l in itertools.product(hO_anos, hi_anos):
                    h1 = myBGF.getAtom(k)   # H atom connected with atom
                    h2 = myBGF.getAtom(l)   # H atom connected with i
                    dist = bgf.distance(h1, h2)
                    if dist > max_dist:
                        max_dist = dist
                        max_dist_pair = [h1, h2]    # now we have two H atoms in maximum distance connected with atom and i

                # dihedral angle
                phi = bgf.dihedral(max_dist_pair[0], atom, i, max_dist_pair[1], radians=False)
                diheds.append(phi)

        hist_angle, _ = np.histogram(angles, bins=bin, normed=True)
        hist_dihed, _ = np.histogram(diheds, bins=bin, normed=True)
        avg_angles.append(hist_angle)
        avg_diheds.append(hist_dihed)

        sys.stdout.write('Done                                                        ')
        sys.stdout.flush()

        # write output
        processed_step += 1
        t2 = time.time()    # time mark
        elapsed_time = t2 - t1;


    avg_angles = np.mean(np.array(avg_angles), axis=0)
    avg_diheds = np.mean(np.array(avg_diheds), axis=0)

    #print(avg_diheds)

    myDAT.write("#angles population\n")
    for index, i in enumerate(avg_angles):
        myDAT.write("%8.2f %8.5f\n" % (_[index], avg_angles[index]))
    myDAT.write("\n\n")
    myDAT.write("#diheds population\n")
    for index, i in enumerate(avg_diheds):
        myDAT.write("%8.2f %8.5f\n" % (_[index], avg_diheds[index]))


    myDAT.close()
    print('')

    return 1

    ### end of function


if __name__ == "__main__":
    bgf_file = ""; trj_file = ""; ff_file = ""; distance = 3.5; condition = ""; nsample = '';
    usage = """
    %s -b bgf_file -t trj_file (-d distance) (-c selection)
        -d distance criteria for the first coordination shell of the water
        -c condition atom selection condition (e.g. atom.chain == "I")
    """ % sys.argv[0]

    options, args = getopt.getopt(sys.argv[1:], 'hb:t:d:c:n:', ['help', 'bgf=', 'trj=', 'distance=', 'condition=', 'nsample='])

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
        elif option in ('-d', '--distance'):
            distance = float(value)
        elif option in ('-c', '--condition'):
            condition = str(value)
        elif option in ('-n', '--nsample'):
            nsample = int(value)
        elif option == NULL:
            print(usage)
            sys.exit(0)
    
    # more options: resname CNT, resname solvent, save BGF

    # main call
    countWaterCNT(bgf_file, trj_file, nsample = nsample, selection = condition, d_crit = distance)
