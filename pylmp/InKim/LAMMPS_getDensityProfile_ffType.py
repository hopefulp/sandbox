#!/home/noische/python

import sys, re, string, getopt, optparse, math, time, pprint 
from os import popen
import bgf
import bgftools
import numpy as np
import operator
import copy
import nutils as nu
import datatools
import lammpstools as lt
#import lammpstrj as trj
from lammps.trj import *
import dreiding
import pickle

option = ""; args = ""; bgf_file = ""; trj_file = ""; out_file = ""; 
version = "160516"

#-----------------
# calculate dipole moments from the LAMMPS trajectory
# http://code.google.com/p/pdb-tools/source/browse/trunk/pdbTools/pdb_moment.py
# last update: 2013/02/06
#_________________
def get_line(trj_file):
    """
    Open a lazy file stream. This is considered to increase a file reading speed.
    """
    with open(trj_file) as file:
        for i in file:
            yield i


def densityProfile(bgf_file, trj_file, out_file, interval, avg_timestep, direction = 'z', ff_file='', silent=False):

    nu.warn("LAMMPS trajectory with NPT simulations will give you the wrong result.")

    ### init
    timestep = 0; l_timestep = []; line = []; n_header = 0;
    t1 = 0; t2 = 0; # clock
    vector = [0, 0, 1];    # the axis of interest. in the acetone-water case, z direction.
    atominfo = bgftools.atom_mass    # atom data extracted from ff_file
    result = dict();
    axis = 0;        # 1: x-axis, 2: y-axis, 3: z-axis


    ### direction
    if "x" in direction:
        axis = 0;
    elif "y" in direction:
        axis = 1;
    elif "z" in direction:
        axis = 2;
    else:
        nu.die("Error on reading direction.")


    ### open files
    myBGF = bgf.BgfFile(bgf_file)
    myTRJ = open(trj_file)
    myTRJ.seek(0)
    f_out_file = open(out_file + ".dat", 'w')
    f_avg_out_file = open(out_file + ".average.dat", 'w')
    f_pickle = open(out_file + ".pickle", 'w')
    f_avg_pickle = open(out_file + ".average.pickle", 'w')


    ### read residues from bgf_file:: residue is replaced by ffTypes
    residue = set();    # kind of residues in BGF file
    dict_residue = dict();    # stores residue numbers per each residue. RES1: [1, 2, ..], RES2: [6, 7, ..]
    for i in myBGF.a:
        rname = string.strip(i.ffType)
        residue.add(rname)
    output = "";
    for i in residue:
        output += i + " ";

    if not silent: print("Found " + str(len(residue)) + " residues in BGF file: " + str(output))
    residue = list(residue)
    residue.append('TOTAL')

    for index, i in enumerate(residue):
        dict_residue[i] = index;

    n_residue = len(residue)    # number of residues (including total)

    if ff_file:
        bgftools.update_mass(ff_file)

    ### read trajectory file
    # how many steps to go?
    #l_timestep = lt.getTrjInfo(trj_file)
    #mytrj = trj.lammpstrj(trj_file)
    mytrj = Trj(trj_file)
    l_timestep = mytrj.load()
    n_timestep = len(l_timestep)
    if not silent: print("\nThe trajectory contains %d timesteps." % n_timestep)

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
        try:
            chunk = [next(dumpatom) for i in range(natoms+n_header)]
        except StopIteration:
            break;
    
        timestep = int(chunk[1])
        l_timestep.append(timestep)
        natoms = int(chunk[3])
        boxsize = [chunk[5].split(' ')[0], chunk[5].split(' ')[1], chunk[6].split(' ')[0], chunk[6].split(' ')[1], chunk[7].split(' ')[0], chunk[7].split(' ')[1]]; 
        boxsize = [float(i) for i in boxsize]; boxsize = [(boxsize[1] - boxsize[0]), (boxsize[3] - boxsize[2]), (boxsize[5] - boxsize[4])]
        keywords = chunk[8].split('ATOMS ')[1].strip('\n').split(' ')

        ### Show progress
        t1 = time.time();
        remaining_time = elapsed_time * (n_timestep - processed_step)
        sys.stdout.write('\r' + "Reading timestep.. " + str(timestep) + " (Elapsed time for the previous step: " + "{0:4.1f}".format(elapsed_time) + " seconds, Remaining time: " + "{0:4.1f}".format(remaining_time) + " seconds = " + "{0:4.1f} minutes".format(remaining_time/60) + ")")
        sys.stdout.flush()

        processed_step += 1;

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
                
        # apply periodic condition
        myBGF = bgftools.periodicMoleculeSort(myBGF, myBGF.CRYSTX, silent=True)

        ### DONE for updating trj on BGF file ###
        ### Do whatever I want: density profile

        ### upper and lower bound
        min = -3.0
        max = boxsize[axis] + 3.0

        ### make a bin
        bin = np.arange(math.floor(min), math.ceil(max), interval)

        ### for every atoms, get atom coordinate and its type
        positions = [];    # atom positions. [ res1, res2, ..., total ]
        masses = [];    # used for a weight for histogram. [ res1, res2, ..., total ]
        hist = [];
        bin_edges = [];
        for i in range(n_residue):
            positions.append([])
            masses.append([])
            hist.append([])
            bin_edges.append([])

        for atom in myBGF.a:
            coord = [atom.x, atom.y, atom.z];
            # total
            positions[-1].append(coord[axis]);
            masses[-1].append(atominfo[atom.ffType]);
            # residues
            positions[dict_residue[atom.ffType.strip()]].append(coord[axis]);
            masses[dict_residue[atom.ffType.strip()]].append(atominfo[atom.ffType]);


        ### calculate volume
        area = 0;
        if axis == 0:
            area = boxsize[1] * boxsize[2]
        elif axis == 1:
            area = boxsize[0] * boxsize[2]
        elif axis == 2:
            area = boxsize[0] * boxsize[1]
        vol = area * interval


        ### histogram
        for i in range(n_residue):
            hist[i], bin_edges[i] = np.histogram(positions[i], bins=bin, weights=masses[i]);
            hist[i] = hist[i] / 6.022 / vol * 10    # density


        ### store
        temp = dict();
        for i in range(n_residue):
            temp2 = dict();
            temp2['HIST'] = hist[i]
            temp2['BIN_EDGES'] = bin_edges[i]
            temp[residue[i]] = temp2;
        result[timestep] = temp


        ### end of loop: check elapsed time
        t2 = time.time()    # time mark
        elapsed_time = t2 - t1;


    ### calculate average values
    if not silent: print("\nAveraging the last %d timesteps.." % avg_timestep)
    avg_timesteps = l_timestep[-avg_timestep:]
    avg_hist = [];
    avg_bin_edges = [];
    for i in range(n_residue):
        avg_hist.append(np.zeros(len(result[l_timestep[-avg_timestep]]['TOTAL']['HIST'])));
        #avg_bin_edges.append(np.zeros(len(result[-avg_timestep]['TOTAL']['BIN_EDGES'])));

    for t in avg_timesteps:
        for r in result[t]:
            avg_hist[dict_residue[r]] += result[t][r]['HIST'];

    for i in range(n_residue):
        avg_hist[i] /= len(avg_timesteps)    # normalize

    temp_avg = dict();    # residue --- HIST, BIN_EDGES
    for i in residue:
        temp2 = dict();
        temp2['HIST'] = avg_hist[dict_residue[i]];
        temp2['BIN_EDGES'] = bin_edges[0]
        temp_avg[i] = temp2


    ### write up a average data to file
    output = str(avg_timestep) + "\n";
    for r in temp_avg:
        output += str(r) + "\n"
        for index, i in enumerate(temp_avg[r]['HIST']):
            output += str(temp_avg[r]['BIN_EDGES'][index]) + "\t" + str(i) + "\n"
        output += "\n"
    f_avg_out_file.write(output)
    

    ### write up a whole data to file
    output = "";
    tkey = result.keys()
    tkey.sort()
    for t in tkey:
        output += str(t) + "\n"
        rkey = result[t].keys()
        rkey.sort()
        for r in rkey:
            output += str(r) + "\n"
            for index, i in enumerate(result[t][r]['HIST']):
                output += str(result[t][r]['BIN_EDGES'][index]) + "\t" + str(i) + "\n"
            output += "\n"
        output += "\n"
    f_out_file.write(output)


    ### write up a pickle object
    if not silent: print("Writing pickle object..")
    pickle.dump(result, f_pickle)
    f_pickle.close()
    pickle.dump(temp, f_avg_pickle)
    f_avg_pickle.close()

    ### return

    print('')
    return 1

    ### end of function


if __name__ == "__main__":
    bgf_file = ""; trj_file = ""; ff_file = ""; out_file = ""; direction = ""; interval = 1; avg_timestep = 0;

    usage = """
    BGF_dipole.py -b bgf_file -t trj_file -f ff_file -o out_file -d direction -i bin_interval -a averaging_timestep
    """

    options, args = getopt.getopt(sys.argv[1:], 'hb:t:f:o:d:i:a:', ['help', 'bgf=', 'trj=', 'ff=', 'out=', 'direction=', 'interval=', 'average='])

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
        elif option in ('-f', '--ff'):
            ff_file = value
        elif option in ('-o', '--out'):
            out_file = value
        elif option in ('-d', '--direction'):
            direction = value
        elif option in ('-i', '--interval'):
            interval = float(value)
        elif option in ('-a', '--average'):
            avg_timestep = int(value)
        elif option == NULL:
            print(usage)
            sys.exit(0)

    direction = string.split(direction)

    if len(direction) == 1:
        if "x" in direction or "y" in direction or "z" in direction:
            pass;
        else:
            nu.die("Please specify the direction to see the density profile.")
    else:
        nu.die("Please specify the correct direction.")

    if out_file == "":
        out_file = bgf_file[:-4] + ".densityProfile";
    
    # main call
    densityProfile(bgf_file, trj_file, out_file, interval, avg_timestep, direction=direction)
