#!/home/noische/python -u

import sys, re, string, getopt, optparse, math, time
from os import popen
import time
import os
import random

import numpy as np

import bgf
import bgftools
import nutils as nu
import lammpstools as lt
import dreiding

option = ""; args = ""; bgf_file = ""; trj_file = ""; out_file = ""; 
usage = """
CNT_traceRadialWaterPosition.py -b bgffile -t trjfile -n #_sample_step -w #_sample_water
"""
version = "160216"

#-----------------
# Get averaged velocity of water molecules confined in CNT
# 
# 
#_________________
def get_line(trj_file):
    """
    Open a lazy file stream. This is considered to increase a file reading speed.
    """
    with open(trj_file) as file:
        for i in file:
            yield i


def countWater(bgf_file, trj_file, n_step, n_water, silent=False):

    ### const
    PI = math.pi
    vdw_r_C = 1.7

    ### init
    timestep = 0; l_timestep = []; line = []; n_header = 0;
    t1 = 0; t2 = 0; # clock

    myBGF = bgf.BgfFile(bgf_file)
    myTRJ = open(trj_file)
    myTRJ.seek(0)

    global out_file
    if out_file == "": out_file = "radialWaterPosition.profile"
    print("The result will be recorded to the file " + out_file + " ...")
    ftemp = open(out_file, 'w')
    ftemp.write(str(sys.argv) + "\n")
    ftemp.write(str(os.path.abspath('.') + "\n"))

    curr_dir = os.path.abspath(".")


    ### how many steps to go?
    n_timestep = len(lt.getTrjInfo(trj_file))
    if n_step == 0:
        n_step = n_timestep;

    print(" ..The trajectory contains " + str(n_timestep) + " timesteps.")
    print("The script will proceed for the first " + str(n_step) + " timesteps.")


    ### extract aNos of CNT in the BGF file
    aNo_CNT = []; aNo_WAT_O = []; 

    for atom in myBGF.a:
        # Carbons in CNT or atoms in BNNT
        if "NT" in atom.rName:
            aNo_CNT.append(atom.aNo)
        # Oxygen in water
        if "WAT" in atom.rName and "O" in atom.aName:
            aNo_WAT_O.append(atom.aNo)

    N_CNT = len(aNo_CNT)    # the number of CNT atoms


    ### Draw O atoms to trace
    n_trace = n_water
    if len(aNo_WAT_O) < n_water:
        nu.die('Cannot draw %d molecules out of %d water molecules in the system!' % (n_water, len(aNo_WAT_O)))
    aNo_WAT_trace = random.sample(aNo_WAT_O, n_trace)
    # convert aNo as rNo
    _ = [];
    for i in aNo_WAT_trace:
        r = myBGF.getAtom(i)
        _.append(r.rNo)
    output = "resid "
    for i in _:
        output += '\t' + str(i)
    output += '\n'
    ftemp.write(output)
    print(output)


    ### check if there exists water properly
    if len(aNo_WAT_O) == 0:
        nu.die("No water molecules in the BGF file.")
    if len(aNo_CNT) == 0:
        nu.die("No CNT molecules in the BGF file.")


    ### Find header of the trajectory file
    while 1:
        templine = myTRJ.readline()
        line.append(templine.strip('\n').strip('ITEM: '))
        n_header += 1
        if "ITEM: ATOMS" in templine:
            break;


    ### INITIAL trajectory information
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
        remaining_time = elapsed_time * (n_step - processed_step)
        sys.stdout.write('\r' + "Reading timestep.. " + str(timestep) + " (Remaining time: " + "{0:4.1f}".format(remaining_time) + " seconds, " + "{0:4.1f} minutes".format(remaining_time/60) + ")")
        sys.stdout.flush()

        if processed_step == n_step:
            break;

        ### Read
        #myBGF = bgf.BgfFile(bgf_file)
        coord = dict()
        try:
            chunk = [next(dumpatom) for i in range(natoms+n_header)]
        except StopIteration:
            break;

        timestep = int(chunk[1])
        natoms = int(chunk[3])
        boxsize = [chunk[5].split(' ')[0], chunk[5].split(' ')[1], chunk[6].split(' ')[0], chunk[6].split(' ')[1], chunk[7].split(' ')[0], chunk[7].split(' ')[1]]; 
        boxsize = [float(i) for i in boxsize]; boxsize = [(boxsize[1] - boxsize[0]), (boxsize[3] - boxsize[2]), (boxsize[5] - boxsize[4])]
        keywords = chunk[8].split('ATOMS ')[1].strip('\n').split(' ')

        mode = 'unwrapped'    # assume that "dump            1 all custom 100 ${sname}${rtemp}K.nvt.lammps id type xu yu zu vx vy vz" in lammps input
    
        # actual coordinate
        coordinfo = chunk[9:]
    
        # load atom coordinates from chunk
        for atomline in coordinfo:
            atomcoord = atomline.split(' ')
            #atom = myBGF.getAtom(int(atomcoord[0]))
            coord[int(atomcoord[0])] = [float(atomcoord[2]), float(atomcoord[3]), 0]
            #atom.x = float(atomcoord[2])
            #atom.y = float(atomcoord[3])
            #atom.z = float(atomcoord[4])

            
        ### convert atom coordinates to lists
        CNT = []; WATER = []; 
        for atom in myBGF.a:
            if "NT" in atom.rName:
                #CNT.append([atom.x, atom.y, atom.z])
                CNT.append(coord[atom.aNo])
            if atom.aNo in aNo_WAT_trace:
                #WATER.append([atom.x, atom.y])
                WATER.append(coord[atom.aNo])

        CNT = np.array(CNT); n_CNT = len(CNT)    # CNT coordinates
        WATER = np.array(WATER)    # chosen water coordinates
        boxsize = np.array(boxsize)

        ### CNT radius
        x_CNT, y_CNT, z_CNT = np.mean(CNT, axis=0)
        CNT_center = np.array([x_CNT, y_CNT, 0])


        ### get radial distance of water
        dist = np.sqrt(((WATER - CNT_center)**2).sum(axis=-1))


        ### write output
        output = str(timestep)
        for i in dist:
            output += '\t' + str(i)
        output += '\n'
        ftemp.write(output)
        sys.stdout.flush()


        ### time
        t2 = time.time()    # time mark
        elapsed_time = t2 - t1;
        processed_step += 1;


    print('')
    ftemp.close()
    print("Numbers of water molecules are written in " + out_file + " ..Done.")


    return 1

    ### end of function


if __name__ == "__main__":
    bgf_file = ""; trj_file = ""; n_step = 0; n_water = 0;

    options, args = getopt.getopt(sys.argv[1:], 'hb:t:n:o:w:', ['help', 'bgf=', 'trj=', 'step=', 'out=', 'nwater='])

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
        elif option in ('-n', '--step'):
            n_step = int(value)
        elif option in ('-o', '--out'):
            out_file = str(value)
        elif option in ('-w', '--nwater'):
            n_water = int(value)
        elif option == NULL:
            print(usage)
            sys.exit(0)
    
    # more options: resname CNT, resname solvent, save BGF

    # main call
    countWater(bgf_file, trj_file, n_step, n_water)
