#!/home/noische/python

import sys
import re
import string
import getopt
import optparse
import time
import os

# BGF module
import bgf
import bgftools
import nutils as nu

# Pizza.py module
#sys.path.append("/home/noische/program/pizza-1Oct10")
#import dump

option = ""; args = ""; 
usage = """
getLAMMPSTrajectory.py: 

Usage: getLAMMPSTrajectory.py -b original_bgf -t trj_file -o out_file -n timestep -p (periodic sort)

"""

version = "160425"

def get_line(trj_file):
    with open(trj_file) as file:
        for i in file:
            yield i

def getLAMMPSTrajectory(bgf_file, trj_file, out_file, request_timestep, periodic, ff_file='', silent=False):
    if os.path.getsize(trj_file) == 0:
        nu.die("Trajectory file size 0. Check trajectory file.")

    if not os.path.isfile(trj_file):
        nu.die("script cannot continue. Check trajectory filename.")

    timestep = 0;
    natoms = 0;
    boxsize = [0, 0, 0, 0, 0, 0];
    keywords = [];
    data = [];
    line = []
    l_timestep_seek = []
    n_header = 0;

    myTrj = open(trj_file)
    myTrj.seek(0)

    # Find and parse header
    while 1:
        templine = myTrj.readline()
        line.append(templine.strip('\n').strip('ITEM: '))
        n_header += 1
        if "ITEM: ATOMS" in templine:
            break;
        
    """
    ### DO NOT DELETE THIS COMMENT ###
    ### LAMMPS Trajectory file header is as follows:

    1: ITEM: TIMESTEP
    2: 0
    3: ITEM: NUMBER OF ATOMS
    4: 5691
    5: ITEM: BOX BOUNDS pp pp pp
    6: 0 44.1947
    7: 0 33.5174
    8: 0 36.6752
    9: ITEM: ATOMS id type x y z ix iy iz
    10: 404 7 6.7879 0.35342 2.62175 0 0 0
    """

    # INITIAL trajectory information to count chunk
    timestep = int(line[1])
    natoms = int(line[3])
    boxsize = [line[5].split(' ')[0], line[5].split(' ')[1], line[6].split(' ')[0], line[6].split(' ')[1], line[7].split(' ')[0], line[7].split(' ')[1]]
    boxsize = [float(i) for i in boxsize]
    keywords = line[8].strip('ATOMS ')

    # Timestep Reading
    dumpatom = get_line(trj_file)
    while 1:
        try:
            chunk = [next(dumpatom) for i in range(natoms+n_header)]

        except StopIteration:
            print("")
            # timestep is not requested from parameters
            if request_timestep != -1:
                nu.warn("Script failed to find the requested timestep.")
            if not silent:    
                print("The last snapshot (timestep " + chunk[1].strip("\n") + ") of the trajectory " + trj_file + " is loaded.")
            break;

        except KeyboardInterrupt:
            print("\nKeyboard Break - Force quit.")
            sys.exit(0)

        else:
            try:
                timestep = int(chunk[1])
            except ValueError:
                print("There was a problem on reading the timestep " + str(timestep))
                print(chunk[1])
                nu.die("Quitting..")
            else:
                if not silent:
                    sys.stdout.write("\r" + "Reading timestep " + chunk[1].strip("\n"))
                    sys.stdout.flush()
                if request_timestep == int(chunk[1]):
                    print("\nRequested timestep found. Timestep loaded.")
                    break;

    # REMARK: now 'chunk' contains the last snapshot information
    # get trajectory information for the requested_timestep
    timestep = int(chunk[1])
    natoms = int(chunk[3])
    boxsize = [chunk[5].split(' ')[0], chunk[5].split(' ')[1], chunk[6].split(' ')[0], chunk[6].split(' ')[1], chunk[7].split(' ')[0], chunk[7].split(' ')[1]]; 
    boxsize = [float(i) for i in boxsize]; boxsize = [(boxsize[1] - boxsize[0]), (boxsize[3] - boxsize[2]), (boxsize[5] - boxsize[4])]
    keywords = chunk[8].split('ATOMS ')[1].strip('\n').split(' ')

    # Open BGF
    if isinstance(bgf_file, bgf.BgfFile):
        myBGF = bgf_file
    else:
        if not silent: print("Reading " + bgf_file + " ..")
        myBGF = bgf.BgfFile(bgf_file)

    natom_bgf = len(myBGF.a)    # number of atoms in BGF file

    if not natom_bgf == natoms:
        nu.die("Number of atoms in trajectory file does not match with BGF file.")

    mode = ""
    if 'xs' in keywords or 'ys' in keywords or 'zs' in keywords:
        if not silent: print("Script found scaled coordinates")
        mode = 'scaled'
    elif 'x' in keywords or 'y' in keywords or 'z' in keywords:
        if not silent: print("Script found normal coordinates")
        mode = 'normal'
    elif 'xu' in keywords or 'yu' in keywords or 'zu' in keywords:
        if not silent: print("Script found unwrapped coordinates")
        mode = 'unwrapped'
    else:
        nu.die("Cannot determine the trajectory dump style.")

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
        if not silent: print("Boxsize: " + str(boxsize))
    except:
        pass;
        

    # apply periodic condition if periodic == True
    if periodic:
        if not silent: print("Applying periodic boundary condition..")
        myBGF = bgftools.periodicMoleculeSort(myBGF, myBGF.CRYSTX, ff_file=ff_file, silent=True, out_file=out_file)

    myBGF.REMARK.insert(0, "Trajectory updated with " + os.path.basename(sys.argv[0]) + " by " + os.environ["USER"] + " on " + time.asctime(time.gmtime()))
    myBGF.REMARK.insert(0, "Trajectory updated with " + str(trj_file))
    myBGF.REMARK.insert(0, "Updated timestep: " + str(timestep))

    # save
    if isinstance(out_file, str):
        if not silent: print("Saving information to " + out_file + " ..")
        myBGF.saveBGF(out_file)
        return 1;
    else:
        return myBGF;

    ## end of getLAMMPSTrajectory


if __name__ == "__main__":

    print("\n" + sys.argv[0] + " version " + str(version) + "\n")

    if len(sys.argv) < 2:
        print(usage);
        sys.exit(1)

    bgf_file = ""; trj_file = ""; out_file = ""; periodic = False; ff_file = ''
    timestep = -1;

    options, args = getopt.getopt(sys.argv[1:], 'hb:t:o:n:f:p', ['help','bgf=','trj=','out=','timestep=','ff=','periodic'])
    try:
        for option, value in options:
            if option in ('-h', '--help'):
                print(usage)
                sys.exit(0);
            elif option in ('-b', '--bgf'):
                bgf_file = value
            elif option in ('-t', '--trj'):
                trj_file = value
            elif option in ('-o', '--out'):
                out_file = value
            elif option in ('-n', '--timestep'):
                timestep = int(value)
            elif option in ('-f', '--ff'):
                ff_file = value
            elif option in ('-p', '--periodic'):
                periodic = True
            elif option == NULL:
                print(usage)
                sys.exit(0)
    except ValueError:
        nu.die("script cannot continue. Check input parameters.")

    if out_file == "":
        out_file = bgf_file.split(".bgf")[0] + ".update.bgf"

    getLAMMPSTrajectory(bgf_file, trj_file, out_file, timestep, periodic, ff_file=ff_file, silent=False)

