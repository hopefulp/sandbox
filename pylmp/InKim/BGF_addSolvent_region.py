#!/home/noische/python

"""
addSolvent.py
Original: Mar 27 2011 In Kim
Version: 160316

# Version updates:

"""

# Python Modules
import sys
import os
import string
import random
import time
import math
import copy

# Custom Modules
sys.path.append("/home/noische/script")
sys.path.append("/home/noische/scripts")
import bgf
import nutils as nu
import tqdm
#from dump import *
import bgftools
from removeBadContacts import *

# Globals
version = '120924'

def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def addsolvent(bgf_file, solvent_bgf, min, max, n_solvent, out_file, ff_file, margin, mark, silent=True):

    ### initialize
    water = False;
    default_margin = 1.0; x_margin = 0.0; y_margin = 0.0; z_margin = 0.0;
    if "x" in margin.lower():
        x_margin = default_margin
    if "y" in margin.lower():
        y_margin = default_margin
    if "z" in margin.lower():
        z_margin = default_margin

    ### load the solute bgf file
    if not silent: print("Initializing..")
    myBGF = bgf.BgfFile(bgf_file)    # str
    myBGF.renumber()

    if not silent: print("Loading the solvent file " + solvent_bgf + " ..")
    solventBGF = bgf.BgfFile(solvent_bgf)

    ### Generate error when the solvent box is not periodic:
    if not solventBGF.PERIOD:
        nu.die("addSolvent: The solvent file is not periodic. Use a box full of solvent.")

    ### Check the type of solvent
    if not silent: print("(the solvent box seems to be full of " + os.path.basename(solvent_file)[:-4] + " )")
    if "f3c" in solvent_bgf: water = True        # this flag is used to remove the bad contacts with the molecule.
    if "spc" in solvent_bgf: water = True        # this flag is used to remove the bad contacts with the molecule.
    if "tip" in solvent_bgf: water = True        # this flag is used to remove the bad contacts with the molecule.

    ### calculate the box size
    if not silent: print("Analyzing box information..")
    if not min:
        min = [0.0, 0.0, 0.0]
    if not max:
        max = myBGF.CRYSTX[:3]

    strsize = [ max[0] - min[0] - x_margin, max[1] - min[1] - y_margin, max[2] - min[2] - z_margin ]

    waterboxsize = solventBGF.CRYSTX[:3]    # REMARK: This is a kind of constant.
    copyNumber = [0, 0, 0];
    for index, i in enumerate(copyNumber):
        copyNumber[index] = math.ceil(strsize[index] / waterboxsize[index])    # how many times to replicate the water box

    if not silent: print("Creating box information: " + str(strsize))

    ### replicate the solvent box
    if not silent: print("Creating box.. this may take some time.")
    bigboxBGF = bgftools.replicateCell(solventBGF, copyNumber, True)
    bigboxBGF.saveBGF("_replicate.bgf")
    if not silent: print("- Number of atoms in the created box: " + str(len(bigboxBGF.a)))
    delatom = []; delsolvent = []; delsolventindex = [];
    delwater = []; delwaterindex = [];

    bigboxBGF = bgf.BgfFile("_replicate.bgf")
    bigboxBGF.renumber()

    ### trim the water box
    if water:
        if not silent: print("Generating water box.. Calculating water molecules")
        for atom in bigboxBGF.a:
            if atom.x > strsize[0] or atom.y > strsize[1] or atom.z > strsize[2]:
                delatom.append(atom.aNo)
        for aNo in delatom:
            water_molecule = bgftools.is_water(bigboxBGF, aNo)
            if not water_molecule in delwater:
                delwater.append(water_molecule)
        delwater = nu.flatten(delwater)
        delwater = nu.removeRepeat(delwater)
        for aNo in delwater:
            delwaterindex.append(bigboxBGF.a2i[aNo])
        if not silent: print("Generating water box.. Trimming")
        bigboxBGF.delAtoms(delwaterindex, False)
        bigboxBGF.renumber()
    elif not water:
        if not silent: print("Generating solvent box.. Extracting solvent molecules")
        for atom in tqdm.tqdm(bigboxBGF.a, desc='Iterating', ncols=120):
            if atom.x > strsize[0] or atom.y > strsize[1] or atom.z > strsize[2]:
                delatom.append(atom.aNo)

        molecule = bgftools.getMoleculeList(bigboxBGF)
        for aNo in tqdm.tqdm(delatom, desc='Appending atoms to remove', ncols=120):
            for i in molecule:
                if aNo in i: 
                    delsolvent += i
                    break
        delsolvent = list(set(delsolvent))
        for aNo in delsolvent:
            delsolventindex.append(bigboxBGF.a2i[aNo])
        if not silent: print("Generating solvent box.. Trimming")
        delsolventindex.sort()
        #delsolventindex.reverse()
        bigboxBGF.delAtoms(delsolventindex, False)
        bigboxBGF.renumber()

    bigboxBGF.saveBGF("_temp.bgf") ## debug

    ### remain n molecules and delete residues
    bigboxBGF = bgftools.renumberMolecules(bigboxBGF, 0, False) # renumber rNo

    if n_solvent:
        residues = set()
        for atom in bigboxBGF.a:
            residues.add(atom.rNo)  # scan molecule rNo
        residues = list(residues)
        if not silent: print("Found %d water molecules." % len(residues))
        if len(residues) < n_solvent:
            nu.die("Too few solvent molecules to choose %d molecules." % n_solvent)
        
        rNos = random.sample(residues, n_solvent)   # select n molecules
        if not silent: print("%d water molecules are chosen." % n_solvent)

        # delete molecules
        delist = []
        for atom in bigboxBGF.a:
            if not atom.rNo in rNos:
                delist.append(bigboxBGF.a2i[atom.aNo])  # if not chosen, delete
        delist.sort()
        bigboxBGF.delAtoms(delist, False);
        bigboxBGF.renumber()

        bigboxBGF = bgftools.renumberMolecules(bigboxBGF, 0, False)

    ### move the water box to min position
    for atom in bigboxBGF.a:
        atom.x += min[0] + x_margin/2
        atom.y += min[1] + y_margin/2
        atom.z += min[2] + z_margin/2

    ### add mark to the solvents
    if mark:
        for atom in bigboxBGF.a:
            atom.chain = mark

    bigboxBGF.saveBGF('_temp.bgf')

    # REMARK: it is natural to have the periodic information of water box for the output BGF file.
    # REMARK: HETATOM should be located on the first of the BGF file. So use dummy for merging.
    ## compute stats for adding solvents
    if not silent: print("\nComputing stats..")
    mol_list = bgftools.getMoleculeList(bigboxBGF)
    n_mol = len(mol_list)
    n_atom = len(nu.flatten(mol_list))
    if not silent: print(str(n_mol) + " molecules (" + str(n_atom) + " atoms) will be added.")

    ## merge
    #bigboxBGF = myBGF.merge(bigboxBGF, True)
    myBGF = bgf.BgfFile(bgf_file)
    bigboxBGF = bgf.BgfFile("_temp.bgf")
    bigboxBGF2 = myBGF.merge(bigboxBGF)
    if not silent: print("Total atoms in the file: " + str(len(bigboxBGF.a)))

    ## some paperworking for periodic box
    bigboxBGF2.OTHER = myBGF.OTHER
    bigboxBGF2.PERIOD = myBGF.PERIOD
    bigboxBGF2.AXES = myBGF.AXES
    bigboxBGF2.SGNAME = myBGF.SGNAME
    bigboxBGF2.CRYSTX = myBGF.CRYSTX
    bigboxBGF2.CELLS = myBGF.CELLS

    ### record BGF remarks
    bigboxBGF2.REMARK.insert(0, "Solvents added by " + os.path.basename(sys.argv[0]) + " by " + os.environ["USER"] + " on " + time.asctime(time.gmtime()))
    bigboxBGF2.REMARK.insert(0, "Solvents: " + str(solvent_file))

    ### save BGF
    if not silent: print("Saving the file.. see " + str(out_file))
    bigboxBGF2.saveBGF(out_file)

    return 1;

    ### end of addSolvent.py


if __name__ == '__main__':
    option = ""; args = ""; bgf_file = ""; ff_file = ""; out_file = ""; solvent_file = ""; ff_file = "";
    suffix = ""; min = ""; max = ""; n_solvent = 0; marginaxis = ''; mark = ''
    usage = """
Usage: addSolvent.py -b bgfFile -f ff_file -m "minx miny minz" -M "maxx maxy maxz" -o outFile -s solventFile
Make the given structure solvated on water box in the given region.

Options:
    -h        print this help message
    -b bgfFile    REQUIRED. The original structure that to be solvated. Only supports BGF format.
    -f ffFile     REQUIRED. CERIUS2 Type Force Field file.
    -o outFile    REQUIRED. Savename.
    -m "n n n"    optional. Minimum x y z points. Same as "0 0 0" (origin) if not specified.
    -M "n n n"    optional. Maximum x y z points. Same as the maximum box boundary if not specified.
    -t solventFile    optional. BGF box filled with solvent. Default is equilibrated F3C water box by Tod Pascal.
            Removing bad contacts btwn HETATM and the solvent will only be supported for F3C water box.
    -n n_solvent    optional. The script will remain n molecules in the system.
    -c "xyz"   optional. Apply 2.0 A margin to the specified axis.
    -v "X"     optional. Change the chain of appending solvents to X.
    """

    if len(sys.argv) < 2:
        print(usage); sys.exit(0)

    options, args = getopt.getopt(sys.argv[1:], 'hb:m:M:o:t:f:n:c:v:', ['help','bgf=','min=','max=','output=','solvent=','ff=','n=','margin=','mark='])
    for option, value in options:
        if option in ('-h', '--help'):
            print(usage); sys.exit(0)
        elif option in ('-b', '--bgf'):
            bgf_file = value
        elif option in ('-m', '--min'):
            min = value
        elif option in ('-M', '--max'):
            max = value
        elif option in ('-o', '--out'):
            out_file = value
        elif option in ('-t', '--solvent'):
            solvent_file = value
        elif option in ('-f', '--ff'):
            ff_file = value
        elif option in ('-n', '--n'):
            n_solvent = int(value)
        elif option in ('-c', '--margin'):
            marginaxis = value
        elif option in ('-v', '--mark'):
            mark = value
        elif option in (''):
            print(usage); sys.exit(0)

    # setting up defaults
    if out_file == "": out_file = bgf_file[:-4] + "_solv" + bgf_file[-4:]
    if solvent_file == "": solvent_file = "/home/noische/scripts/dat/WAT/spc_box.bgf"
    if ff_file == "": nu.die("No Force Field File!")

    # min & max defaults
    mybgf = bgf.BgfFile(bgf_file)
    pbc = mybgf.CRYSTX[:3]

    if min:
        value = min.split()
        if len(value) != 3:
            nu.die("Lower bound not properly set.")

        min = []
        for i in value:
            if is_float(i):
                min.append(float(i))
            else:
                min.append(0)
    else:
        min = []  # no min specified

    if max:
        value = max.split()
        if len(value) != 3:
            nu.die("Upper bound not properly set.")

        max = []
        for index, i in enumerate(value):
            if is_float(i):
                max.append(float(i))
            else:
                max.append(pbc[index])
    else:
        max = []  # no max specified

    # mark should be one letter
    if len(mark) > 2: nu.die("Chain Mark (option -v) should be a character.")

    # main call without silent
    addsolvent(bgf_file, solvent_file, min, max, n_solvent, out_file, ff_file, marginaxis, mark, False)
