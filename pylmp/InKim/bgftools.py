"""
bgftools.py
Original: Jan 01 2011 In Kim

Module containing BGF-file information extraction tools including
getCom(myBGF, forcefield):
centerBGF(myBGF):
"""

# python modules
import sys
import os
import string
import copy
import math
import operator
from types import *

# custom modules
import bgf
import nutils as nu
import dreiding
import numpy as np

version = '111206'

# constants
# from DREIDING2.22.ff
atom_mass = {'H_': 1.00800, \
            'H___A': 1.00800, \
            'H___b': 1.00800, \
            'B_3':  10.81000, \
            'B_2':  10.81000, \
            'C_34': 16.04300, \
            'C_33': 15.03500, \
            'C_32': 14.02700, \
            'C_31': 13.01900, \
            'C_3':  12.01100, \
            'C_2G': 12.01100, \
            'C_22': 14.02700, \
            'C_21': 13.01900, \
            'C_2':  12.01100, \
            'C_R2': 14.02700, \
            'C_R1': 13.01900, \
            'C_R':  12.01100, \
            'C_11': 13.01900, \
            'C_1':  12.01100, \
            'N_3':  14.00670, \
            'N_2':  14.00670, \
            'N_R':  14.00670, \
            'N_1':  14.00670, \
            'O_3':  15.99940, \
            'O_2':  15.99940, \
            'O_R':  15.99940, \
            'O_1':  15.99940, \
            'F_':   18.99840, \
            'Al3':  26.98150, \
            'Si3':  28.08600, \
            'P_3':  30.97380, \
            'S_3':  32.06000, \
            'Cl':   35.45300, \
            'Ga3':  69.72000, \
            'Ge3':  72.59000, \
            'As3':  74.92160, \
            'Se3':  78.96000, \
            'Br':   79.90400, \
            'In3': 114.82000, \
            'Sn3': 118.69000, \
            'Sb3': 121.75000, \
            'Te3': 127.60000, \
            'I_':  126.90450, \
            'Na':   22.99000, \
            'Ca':   40.08000, \
            'Ti':   47.90000, \
            'Fe':   55.84700, \
            'Zn':   65.37700, \
            'Tc':   98.90620, \
            'Ru':  101.07000, \
            'Mo':   95.95000, \
            'S_3a': 32.06000, \
            'S_3b': 32.06000, \
            'O_3F': 15.99940, \
            'H_F':   1.00800, \
            'OW':   16.00000, \
            'HW':    1.00800 }

loaded_ff = ['DREIDING2.21.ff']; # used to memorize already loaded ff_files.
loaded_types = [];

def update_mass(ff_file, silent=True):
    for i in ff_file:
        # no need to read the ff_file again.
        if i in loaded_ff:
            continue;

        # read
        if not silent: print("Updating atom mass with %s" % i)
        FF = dreiding.loadFF(i)
        atominfo = dreiding.loadAtomTypes(FF)
        for key in atominfo.keys():
            atom_mass[key] = atominfo[key]['MASS']
        loaded_ff.append(i)


def get_ffTypes(myBGF, silent=False):
    fftypes = set()
    for atom in myBGF.a:
        fftypes.add(atom.ffType)

    return list(fftypes)


def getCom(myBGF, ff_file = '', aNo_list = [], silent = False):
    """
getCom(myBGF):
    Returns a (x, y, z) of the center of mass of the molecule.
    """

    # init
    mrx = 0; mry = 0; mrz = 0; m = 0;

    # if aNo_list not specified, CoM of whole molecule will be returned.
    if not aNo_list:
        for atom in myBGF.a:
            aNo_list.append(atom.aNo)

    # read atom mass from the dictionary
    for i in aNo_list:
        atom = myBGF.getAtom(i)
        fftype_key = atom.ffType.strip()
        if not fftype_key in atom_mass:
            parse = ff_file.split()
            update_mass(parse, silent=silent)

        if not fftype_key in atom_mass:
            nu.die("No ffType %s found in the atom mass dictionary. You have to specify a ff_file to countinue." % fftype_key)

        mrx += ( atom.x * float(atom_mass[fftype_key]) )
        mry += ( atom.y * float(atom_mass[fftype_key]) )
        mrz += ( atom.z * float(atom_mass[fftype_key]) )
        m += float(atom_mass[fftype_key])

    try:
        x = mrx / m; y = mry / m; z = mrz / m;
    except:
        nu.warn("Total mass is zero for the molecule while calculating center of mass!")
        x = 0; y = 0; z = 0;

    return (x, y, z)


def getMass(myBGF, aNo_list = [], ff_file = "", silent = False):
    """
def getMass(myBGF, aNo_list, ff_file):
    Returns a sum of atom masses in aNo_list.
    """

    if ff_file:
        parse = ff_file.split()
        update_mass(parse)

    total_mass = 0;

    if not aNo_list:
        for atom in myBGF.a:
            aNo_list.append(atom.aNo)

    for ano in aNo_list:
        atom = myBGF.getAtom(ano)
        fftype = atom.ffType.strip()
        try:
            total_mass += float(atom_mass[fftype])
        except KeyError:
            nu.warn("No atom mass found for atomtype %s!" % fftype)
            total_mass += 0.0;

    return total_mass;


possible_atom_attr = ['atom.x', 'atom.y', 'atom.z', 'atom.charge', 'atom.vx', 'atom.vy', 'atom.vz', 'atom.fx', 'atom.fy', 'atom.fz']
def atoms_average(mybgf, attr, selection=''):
    '''
    attr: e.g.) x, y, z or atom.x, atom.y, ...
    selection: e.g.) "'Mo' in atom.ffType"
    '''
    # init
    avg = 0.0; n = 0;

    # conditions
    if not attr:
        nu.die("You have to specify BgfAtom attribute to average!")

    if not 'atom.' in attr:
        attr = "atom." + attr

    if not attr in possible_atom_attr:
        nu.die("Impossible to get average of %s!" % attr)

    if not "atom" in selection:
        nu.die("You have to specify a proper selection!")

    # get average
    for atom in mybgf.a:
        if selection:
            if eval(selection):
                avg += eval(attr)
                n += 1
        else:
            avg += eval(attr)
            n += 1

    return avg/float(n)


def atoms_minmax(mybgf, attr, selection=''):
    min = 1e10; max = -1e10;

    # conditions
    if not attr:
        nu.die("You have to specify BgfAtom attribute to get minmax!")

    if not 'atom.' in attr:
        attr = "atom." + attr

    if not attr in possible_atom_attr:
        nu.die("Impossible to get minmax of %s!" % attr)

    if not "atom" in selection:
        nu.die("You have to specify a proper selection!")

    # get minmax
    for atom in mybgf.a:
        if selection:
            if eval(selection):
                if eval(attr) < min:
                    min = eval(attr)
                if eval(attr) > max:
                    max = eval(attr)
        else:
            if eval(attr) < min:
                min = eval(attr)
            if eval(attr) > max:
                max = eval(attr)

    return min, max


def get_neighbors_aNo(atoms, center, r, pbc="", k=6, return_distance=False):
    """
    returns a list of atom numbers within r of the point.
    atoms: [bgf.BgfAtoms]
    center: np.array
    r = float
    """
    coords = []; aNos = []; neighbors = [];

    for atom in atoms:
        coords.append([atom.x, atom.y, atom.z])
        aNos.append(atom.aNo)

    coords = np.array(coords)

    if pbc:
        import pkdtree
        tree = pkdtree.PeriodicKDTree(pbc, coords, leafsize=len(coords)+1)
    else:
        from scipy.spatial import KDTree
        tree = KDTree(coords, leafsize=len(coords)+1)

    dist, ndx = tree.query([center], k=k)

    if center in coords:
        index = np.where((dist >= 1e-6) & (dist <= r)) # discard center itself
    else:
        index = np.where(dist <= r)

    for i in ndx[index]:
        neighbors.append(aNos[i])

    if not return_distance:
        return neighbors
    else:
        return neighbors, dist


def is_hbonded(mybgf, atom1, atom2, d_crit=3.5, a_crit=30.0):
    """
    returns True if two atoms are h-bonded with a popular geometric definition d(O..O) = 3.5 and A(H-O-O) < 30'
    20160529: Works only for O atoms.
    """
    from numpy import arccos
    from numpy.linalg import norm

    if not "O" in atom1.ffType or not"O" in atom2.ffType:
        return False

    d = np.array([atom1.x, atom1.y, atom1.z])
    a = np.array([atom2.x, atom2.y, atom2.z])
    if 1e-5 < nu.pbc_dist(a, d, mybgf.CRYSTX[:3]) < d_crit:
        for ano in atom1.CONECT:
            # atom1: donor, atom2: acceptor
            h_atom = mybgf.getAtom(ano)
            h = np.array([h_atom.x, h_atom.y, h_atom.z])
            u = h - d; v = a - d;
            theta = np.dot(u, v) / norm(u) / norm(v); theta = np.degrees(arccos(theta))
            if theta < a_crit:
                return True

        for ano in atom2.CONECT:
            # atom2: donor, atom1: acceptor
            h_atom = mybgf.getAtom(ano)
            h = np.array([h_atom.x, h_atom.y, h_atom.z])
            u = h - d; v = a - d;
            theta = np.dot(u, v) / norm(u) / norm(v); theta = np.degrees(arccos(theta))
            if theta < a_crit:
                return True

    return False


def is_hbonded2(mybgf, atom1, atom2, d_crit=2.2):
    """
    returns True if two atoms are h-bonded by a geometric definition.
    This definition is introduced in Grishina, N., & Buch, V. (2004). J. Chem. Phys., 120(11), 5217.
    20160608: Works only for O atoms.
    """
    
    if not "O" in atom1.ffType or not"O" in atom2.ffType:
        return False

    d = np.array([atom1.x, atom1.y, atom1.z])
    a = np.array([atom2.x, atom2.y, atom2.z])

    if mybgf.CRYSTX:
        pbc = mybgf.CRYSTX[:3]
    else:
        pbc = 0

    # A1-H ... A2
    for ano in atom1.CONECT:
        atom = mybgf.getAtom(ano)   # this should be an H atom
        x = np.array([atom.x, atom.y, atom.z])
        if pbc:
            dist = nu.pbc_dist(x, a, pbc)
        else:
            dist = nu.dist(x, a)
        if 1e-5 < dist < d_crit:
            return x

    # A1 ... H-A2
    for ano in atom2.CONECT:
        atom = mybgf.getAtom(ano)   # this should be an H atom
        x = np.array([atom.x, atom.y, atom.z])
        if pbc:
            dist = nu.pbc_dist(x, a, pbc)
        else:
            dist = nu.dist(x, a)
        if 1e-5 < dist < d_crit:
            return x

    return []


def is_water(myBGF, index):
    """
is_water(myBGF, index):
    check if this atom is contained in water molecule or not and returns the aNo of a water molecule.
    """
    dest_atom = myBGF.getAtom(index);

    # for Oxygen
    if "O" in dest_atom.ffType:
        # it has only two atoms
        if len(dest_atom.CONECT) == 2:
            # ..and they are all hydrogens
            if "H" in myBGF.getAtom(dest_atom.CONECT[0]).ffType and "H" in myBGF.getAtom(dest_atom.CONECT[1]).ffType:
                # ..and their hydrogens has no connected atoms
                #if len(myBGF.getAtom(dest_atom.CONECT[0]).CONECT) == 1 and len(myBGF.getAtom(dest_atom.CONECT[1]).CONECT) == 1:
                #    return [dest_atom.aNo, myBGF.getAtom(dest_atom.CONECT[0]).aNo, myBGF.getAtom(dest_atom.CONECT[1]).aNo]
                return [dest_atom.aNo, myBGF.getAtom(dest_atom.CONECT[0]).aNo, myBGF.getAtom(dest_atom.CONECT[1]).aNo]
    # for Hydrogen
    elif "H" in dest_atom.ffType:
        alist = dest_atom.CONECT
        if len(alist) == 1:
            if is_water(myBGF, myBGF.getAtom(alist[0]).aNo):
                # O H1 H2
                return [myBGF.getAtom(alist[0]).aNo, myBGF.getAtom(alist[0]).CONECT[0], myBGF.getAtom(alist[0]).CONECT[1]]
            else:
                return [];
        else:
            return [];
    else:
        return [];

    return [];


def listAllAtoms(myBGF):
    """
listAllAtoms(myBGF):
    returns a list of all atom index.
    """
    atom_list = []
    for atom in myBGF.a:
        atom_list.append(atom.aNo)

    return atom_list;


def listOxygenAtoms(myBGF):
    """
listOxygenAtoms(self, myBGF):
    returns a list of oxygen atom index.
    """
    oxygen_list = [];
    for atom in myBGF.a:
        if "O" in atom.ffType:
            oxygen_list.append(atom.aNo)

    return oxygen_list;


def listWaterAtoms(myBGF):
    """
listWaterAtoms(self, myBGF):
    returns a whole list of oxygen and hydrogen atom indices which is in water.
    """
    water_list = [];

    for atom_index in listOxygenAtoms(myBGF):
        if is_water(myBGF, atom_index):
            water = myBGF.getAtom(atom_index);
            water_list.append(atom_index)
            for conect in water.CONECT:
                water_list.append(conect)

    return water_list;


def listSoluteAtoms(myBGF):
    """
listSoluteAtoms(myBGF):
    returns a list of solute atoms (i.e. system - water).
    """
    all_list = listAllAtoms(myBGF)

    for water_aNo in listWaterAtoms(myBGF):
        if water_aNo in all_list:
            all_list.remove(water_aNo)

    return all_list;


def deleteWaterAtoms(myBGF, index):
    """
deleteWaterAtoms(myBGF, index):
    deletes a water molecule which contains an aNo index.
    """

    water = is_water(myBGF, index)    # find the water molecule which contains the specified oxygen or hydrogen atom.
    water_aNo = []
    for i in water:
        water_aNo.append(myBGF.getAtomIndex(i))        # make a INDEX list for removal (not aNo)

    #print("Deleting the water molecule: " + str(water_aNo))
    myBGF.delAtoms(water_aNo)        # delete the water molecule according to their INDEX (not aNo)
    myBGF.renumber()

    return myBGF


def replicateCell(myBGF, multiplication, spaceFilling = True):
    """
replicateCell():
    Requires a BgfFile class.
    Returns a BgfFile class of the replicated structure of the given structure.

    Options:
    spaceFilling = True
        Just copy the cell (x, y, z) times to the (x, y, z) coordinateas. Negative multiplication indices are allowed.
        (1, 1, 1) with a cell gives 4 cells (i.e. replication to x, y, and z coordinates).

    spaceFilling = False
        Fill the space as a Cuboid form for a given times. Negative multiplication indices are NOT allowed.
        (2, 1, 3) with a cell gives the 6 cells (i.e. cubiod). (1, 1, 1) with a cell gives the identical structure (nothing happens).
    """

    copyTime = [int(multiplication[0]), int(multiplication[1]), int(multiplication[2])];        # the number of cells that will be copied
    myBGFx = bgf.BgfFile(); myBGFy = bgf.BgfFile(); myBGFz = bgf.BgfFile();
    boxsize = [0, 0, 0];

    if len(myBGF.CRYSTX) > 2:
        boxsize = [ myBGF.CRYSTX[0], myBGF.CRYSTX[1], myBGF.CRYSTX[2] ]
    else:
        box = bgf.getBGFSize(myBGF, 0)          # [xlo, xhi, ylo, yhi, zlo, zhi]
        boxsize = [ box[1] - box[0], box[3] - box[2], box[5] - box[4] ]

    if spaceFilling:
        if copyTime[0] != 0:
            myBGF2 = copy.deepcopy(myBGF)
            for n in range(0, copyTime[0] - 1):
                bgf.moveBGF(myBGF2, boxsize[0], 0, 0)
                myBGF = myBGF.merge(myBGF2, True)
        if copyTime[1] != 0:
            myBGF2 = copy.deepcopy(myBGF)
            for n in range(0, copyTime[1] - 1):
                bgf.moveBGF(myBGF2, 0, boxsize[1], 0)
                myBGF = myBGF.merge(myBGF2, True)
        if copyTime[2] != 0:
            myBGF2 = copy.deepcopy(myBGF)
            for n in range(0, copyTime[2] - 1):
                bgf.moveBGF(myBGF2, 0, 0, boxsize[2])
                myBGF = myBGF.merge(myBGF2, True)

        # update the CRYSTX information
        if len(myBGF.CRYSTX) > 2:
            for index, number in enumerate(copyTime):
                myBGF.CRYSTX[index] *= copyTime[index]

    else:
        # replication
        if copyTime[0] != 0:
            if copyTime[0] < 0: sign = -1
            else: sign = 1
            for n in range(0, abs(copyTime[0])):
                if copyTime[0] < 0: n -= 1
                myBGF2 = copy.deepcopy(myBGF)
                bgf.moveBGF(myBGF2, boxsize[0] * (n + 1) * sign, 0, 0)
                myBGFx = myBGFx.merge(myBGF2, True)
        if copyTime[1] != 0:
            if copyTime[0] < 0: sign = -1
            else: sign = 1
            for n in range(0, abs(copyTime[1])):
                if copyTime[0] < 0: n -= 1
                myBGF2 = copy.deepcopy(myBGF)
                bgf.moveBGF(myBGF2, 0, boxsize[1] * (n + 1) * sign, 0)
                myBGFy = myBGFy.merge(myBGF2, True)
        if copyTime[2] != 0:
            if copyTime[0] < 0: sign = -1
            else: sign = 1
            for n in range(0, abs(copyTime[2])):
                if copyTime[0] < 0: n -= 1
                myBGF2 = copy.deepcopy(myBGF)
                bgf.moveBGF(myBGF2, 0, 0, boxsize[2] * (n + 1) * sign)
                myBGFz = myBGFz.merge(myBGF2, True)
        
        myBGF = myBGF.merge(myBGFx, True)
        myBGF = myBGF.merge(myBGFy, True)
        myBGF = myBGF.merge(myBGFz, True)

        # update the CRYSTX information
        if len(myBGF.CRYSTX) > 2:
            for index, number in enumerate(copyTime):
                myBGF.CRYSTX[index] += myBGF.CRYSTX[index] * copyTime[index]
        
    return myBGF


def getmolecule(myBGF, atom, list_aNo, recursion_limit=10000):
    """
getmolecule(myBGF, atom, list_aNo):
    returns a list of atom numbers (aNo) which is linked to the designated atom number.

    myBGF:    bgf.BgfFile class
    atom:    bgf.BgfAtom class
    list_aNo:    a "BLANK" list (IMPORTANT! this should be declared)
    """

    sys.setrecursionlimit(recursion_limit)

    if atom.aNo not in list_aNo:
        list_aNo.append(atom.aNo)
    else:
        return list_aNo;

    for i in atom.CONECT:
        nextatom = myBGF.getAtom(i)
        getmolecule(myBGF, nextatom, list_aNo)


def getBackbone(bgf_file, out_file=False, silent=True):
    """
def getBackbone(string, string): return 1
def getBackbone(BgfFile): return BgfFile
    returns a dehydrogenated BGF structure.

Function Parameters:
    bgf_file    A string which contains a monomer information in a BGF format. 
    """

    # initialize
    temp_index = [];

    # open bgf
    if isinstance(bgf_file, bgf.BgfFile):
        myBGF = bgf_file
    else:
        if not silent: print("reading " + bgf_file + " ..")
        myBGF = bgf.BgfFile(bgf_file)

    # repeat up to number:
    if not silent: print("deleting all ATOMs and hydrogens..")
    for atom in myBGF.a:
        if atom.aTag == 0:    # if atom is an ATOM
            temp_index.append(myBGF.a2i[atom.aNo])
        if "H_" in atom.ffType:
            temp_index.append(myBGF.a2i[atom.aNo])

    myBGF.delAtoms(temp_index)
    myBGF.removeDanglingBonds()
    myBGF.renumber()
    ### REMARK: renumbering should not be done to sustain appropriate aNos.

    # save
    if isinstance(out_file, str):
        if not silent: print("saving information to " + out_file + " ..")
        myBGF.saveBGF(out_file)
        return 1;
    else:
        return myBGF;

    ### end of backbone


def getConnectionDict(bgf_file, silent=True):
    """
def getConnectionDict(str):
def getConnectionDict(BgfFile):
    Returns a dictionary of atom connection for the input of findShortestPath()

Function Parameters:
    bgf_file:    a filename or a BgfFile class (automatically checked)
    """
    # open bgf
    if isinstance(bgf_file, bgf.BgfFile):
        myBGF = bgf_file
    else:
        if not silent: print("reading " + bgf_file + " ..")
        myBGF = bgf.BgfFile(bgf_file)

    connection = dict();
    for atom in myBGF.a:
        connection[atom.aNo] = atom.CONECT

    return connection;

    ### end of getConnectionDict


def findShortestPath(graph, start, end, path=[]):
    """
def findShortestPath(graph, start, end):
    Using backtracking algorithm: from http://www.python.org/doc/essays/graphs.html

Function Parameters:
    graph: a dictionary created by getConnectionDict()
    start: starting atom number
    end: ending atom number
    """
    # recursion setting
    #sys.setrecursionlimit(10000)

    path = path + [start]
    if start == end:
        return path
    if not graph.has_key(start):
        return None
    shortest = None
    for node in graph[start]:
        if node not in path:
            newpath = findShortestPath(graph, node, end, path)
            if newpath:
                if not shortest or len(newpath) < len(shortest):
                    shortest = newpath
    return shortest

    ### end of findShortestPath


def getShortestPath(bgf_file, ano1, ano2, silent=True):
    """
def getShortestPath(bgf_file, ano1, ano2):
    Get a list of the shortest path between two atoms of the given BGF structure using findShortestPath()

Function Parameters:
    bgf_file:    a filename or BgfFile class
    ano1:        atom number of the starting atom
    ano2:        atom number of the ending atom
    """

    # open bgf
    if isinstance(bgf_file, bgf.BgfFile):
        myBGF = bgf_file
    else:
        if not silent: print("reading " + bgf_file + " ..")
        myBGF = bgf.BgfFile(bgf_file)

    ### REMARK: be sure that the hydrogen branches of bgf_file are taken off.
    ## cut off all hydrogen branches
    #myBGF = getBackbone(myBGF)

    # prepare the connection information
    graph = getConnectionDict(myBGF)

    # run and return
    return findShortestPath(graph, ano1, ano2);

    ### end of getShortestPath


def getMoleculeList(bgf_file, recursion_limit=10000, silent=True):
    """
def getMoleculeList(bgf_file, silent=True):
    Get a list of atom numbers (aNo) of the molecules of the given BGF file.
    """

    sys.setrecursionlimit(recursion_limit)

    # open bgf
    if isinstance(bgf_file, bgf.BgfFile):
        myBGF = bgf_file
    else:
        if not silent: print("reading " + bgf_file + " ..")
        myBGF = bgf.BgfFile(bgf_file)

    # build all ano list
    l_all_atom = [];
    for atom in myBGF.a:
        l_all_atom.append(atom.aNo)

    # for every atoms, get connected atoms
    l_all_molecules = [];
    l_molecule = [];
    while len(l_all_atom) > 0:
        l_temp = []; 
        l_molecule = getmolecule(myBGF, myBGF.getAtom(l_all_atom[0]), l_temp, recursion_limit)    # l_molecule is kind of dummy.. no items in there.
        l_all_molecules.append(l_temp)
        for i in l_temp:
            l_all_atom.remove(i)
    
    return l_all_molecules;

    ### end of getMoleculeList


def renumberMolecules(bgf_file, out_file, silent=True):
    """
def renumberMolecules(bgf_file, out_file, silent=True):
    reassign residue number in the molecular order (e.g. water molecules by 1, 2, 3, ...)
    """

    # open bgf
    if isinstance(bgf_file, bgf.BgfFile):
        myBGF = bgf_file
    else:
        if not silent: print("Reading " + bgf_file + " ..")
        myBGF = bgf.BgfFile(bgf_file)

    # get lists of molecules
    l_molecule = getMoleculeList(bgf_file)
    l_molecule = sorted(l_molecule, key=lambda x:sorted(x)[0])
    nmol = len(l_molecule)
    if not silent: print(str(nmol) + " Fragments exists in the BGF file.")
    if nmol > 1000: print("Warning: more than 1000 atoms can cause unwanted errors while treating the BGF file.")

    molnum = 1;
    # for every molecules
    if not silent: print("Renumbering molecules..")
    for molecule in l_molecule:
        # change their molecules number
        for aNo in molecule:
            myBGF.getAtom(aNo).rNo = molnum
        molnum += 1;

    # save
    if isinstance(out_file, str):
        if not silent: print("Saving information to " + out_file + " ..")
        myBGF.saveBGF(out_file)
        return 1;
    else:
        return myBGF;

    ### end of renumberMolecules


def periodicMoleculeSort(bgf_file, boxinfo, fragments=[], selection='', out_file='', ff_file="", silent=False, recursion_limit=10000):
    """
def periodicMoleculeSort(bgf_file, out_file, boxinfo, silent=True, recursion_limit=10000):

    read a BGF format from bgf_file and write the periodic information to out_file
    fragments = [ [atom no. of molecule 1], [atom no. of molecule 2], ... ]
    """
    # open bgf
    if isinstance(bgf_file, bgf.BgfFile):
        if not silent: print("bgftools.periodicMoleculeSort: Got a tossed BGF format ..")
        myBGF = bgf_file
    else:
        if not silent: print("bgftools.periodicMoleculeSort: Reading " + bgf_file + " ..")
        myBGF = bgf.BgfFile(bgf_file)
    
    # if there is a requesting box size, assign
    if boxinfo:
        myBGF.CRYSTX = boxinfo;    # assign requsted boxinfo to CRYSTX

    # get periodic information of bgf file
    l_myBGFbox = [];    # box size
    if myBGF.CRYSTX == []:
        if not silent: 
            nu.die("bgftools.periodicMoleculeSort: ERROR: This BGF file is not periodic. Exiting.")
        else:
            print("bgftools.periodicMoleculeSort: ERROR: This BGF file is not periodic. Exiting.")
            sys.exit(0)
    else:
        l_myBGFbox = copy.deepcopy(myBGF.CRYSTX)

    # get molecules information
    if fragments:
        l_all_molecules = fragments
    else:
        l_all_molecules = getMoleculeList(myBGF, recursion_limit)

    if selection:
        selected_atoms = [atom.aNo for atom in myBGF.a if eval(selection)]
        new_molecules_list = [mol for mol in l_all_molecules for ano in selected_atoms if ano in mol]
        new_molecules_list = nu.removeRepeat(new_molecules_list)
    else:
        new_molecules_list = l_all_molecules

    # for every molecules..
    for molecule in new_molecules_list:
        # calculate COM
        cx, cy, cz = getCom(myBGF, ff_file=ff_file, aNo_list = molecule, silent=True)

        # if COM is out of the box then take the residue
        qx, _ = divmod(cx, l_myBGFbox[0])
        qy, _ = divmod(cy, l_myBGFbox[1])
        qz, _ = divmod(cz, l_myBGFbox[2])
        for aNo in molecule:
            atom = myBGF.getAtom(aNo)
            atom.x = atom.x - l_myBGFbox[0] * qx
            atom.y = atom.y - l_myBGFbox[1] * qy
            atom.z = atom.z - l_myBGFbox[2] * qz

    # save
    if out_file:
        if not silent: print("saving information to " + out_file + " ..")
        myBGF.saveBGF(out_file)

    return myBGF;

    ### end of periodicMoleculeSort


def getAmineGroupInfo(bgf_file, silent=True):
    """
getAmineGroupInfo(bgf_file, silent=True):
    get a BGF file or BgfFile class object
    returns a triplet list of numbers of primary, secondary, and tertiary amines.

    Note that any corrections are not applied to the original BGF file.
    """
    n_pri = 0; n_sec = 0; n_ter = 0; n_garbage = 0;

    # open bgf
    if isinstance(bgf_file, bgf.BgfFile):
        myBGF = bgf_file
    else:
        if not silent: print("Reading " + bgf_file + " ..")
        myBGF = bgf.BgfFile(bgf_file)

    # get nitrogen ano lists
    l_nitrogen_ano = [];
    for atom in myBGF.a:
        if "N_" in atom.ffType: l_nitrogen_ano.append(atom.aNo)

    # for every nitrogen
    for aNo in l_nitrogen_ano:
        n_carbon = 0; 
        connected_ano = myBGF.getAtom(aNo).CONECT
        ## for every connection
        for aNo2 in connected_ano:
            if "C_" in myBGF.getAtom(aNo2).ffType: 
                n_carbon += 1

        ### no carbon: primary, 1 carbons: secondary, 2 carbons: tertiary
        if n_carbon == 1: 
            n_pri += 1;
            #myBGF.getAtom(aNo).rName = "PRI"
        elif n_carbon == 2:
            n_sec += 1;
            #myBGF.getAtom(aNo).rName = "SEC"
        elif n_carbon == 3:
            n_ter += 1;
            #myBGF.getAtom(aNo).rName = "TER"
        else:
            n_garbage += 1;

    # return the number of pri, sec, and ter
    return [n_pri, n_sec, n_ter]

    ### end of getAmineGroupInfo


def setAmineGroupInfo(bgf_file, out_file, silent=True):
    """
setAmineGroupInfo(bgf_file, silent=True):
    get a BGF file or BgfFile class object
    returns a BgfFile or 1 with  primary, secondary, and tertiary amines information corrected.

    Note that corrections WILL BE applied to the original BGF file.
    """
    n_pri = 0; n_sec = 0; n_ter = 0; n_garbage = 0;

    # open bgf
    if isinstance(bgf_file, bgf.BgfFile):
        myBGF = bgf_file
    else:
        if not silent: print("Reading " + bgf_file + " ..")
        myBGF = bgf.BgfFile(bgf_file)

    # get nitrogen ano lists
    l_nitrogen_ano = [];
    for atom in myBGF.a:
        if "N_" in atom.ffType: l_nitrogen_ano.append(atom.aNo)

    # for every nitrogen
    for aNo in l_nitrogen_ano:
        n_carbon = 0; 
        connected_ano = myBGF.getAtom(aNo).CONECT
        ## for every connection
        for aNo2 in connected_ano:
            if "C_" in myBGF.getAtom(aNo2).ffType: 
                n_carbon += 1

        ### no carbon: primary, 1 carbons: secondary, 2 carbons: tertiary
        if n_carbon == 1: 
            myBGF.getAtom(aNo).rName = "PRI"
        elif n_carbon == 2:
            myBGF.getAtom(aNo).rName = "SEC"
        elif n_carbon == 3:
            myBGF.getAtom(aNo).rName = "TER"
        else:
            n_garbage += 1;    # what are you??

    if n_garbage != 0:
        nu.warn("Suspicious amine groups are found!")
        return 0;

    # save
    if isinstance(out_file, str):
        if not silent: print("Saving information to " + out_file + " ..")
        myBGF.saveBGF(out_file)
        return getAmineGroupInfo(myBGF);
    else:
        return getAmineGroupInfo(myBGF);

    ### end of setAmineGroupInfo


def removeFragments(bgf_file, resname, n_frag, out_file, silent=True):
    """
    """

    # variables
    l_delatoms = []; n_molecules = 0; n_atoms = 0;

    # open bgf
    if isinstance(bgf_file, bgf.BgfFile):
        myBGF = bgf_file
    else:
        if not silent: print("Reading " + bgf_file + " ..")
        myBGF = bgf.BgfFile(bgf_file)

    # get lists of molecules
    l_molecule = getMoleculeList(bgf_file)

    if not silent: print("Molecules which is less atoms than " + str(n_frag) + " among the residue name " + str(resname) + " will be removed.")
    
    # for every molecule
    for molecule in l_molecule:
        if len(molecule) < n_frag:
            for ano in molecule: 
                # check rName
                atom = myBGF.getAtom(ano)
                if resname in atom.rName:
                    l_delatoms.append(myBGF.a2i[ano])    # remove index, not ano!
                    #l_delatoms.append(ano)    # remove index, not ano!
                    n_molecules += 1

    # delete and compute stats
    n_atoms = len(l_delatoms)
    l_delatoms.sort()
    l_delatoms.reverse()
    myBGF.renumber()
    myBGF.delAtoms(l_delatoms)
    myBGF.renumber()

    if not silent: print(str(n_molecules) + " molecules (" + str(n_atoms) + " atoms) are deleted.")

    # save
    if isinstance(out_file, str):
        if not silent: print("Saving information to " + out_file + " ..")
        myBGF.saveBGF(out_file)
        return 1;
    else:
        return myBGF;

    ### end of removeFragments


def selectAtoms(bgf_file, selection, out_file, silent=True):
    """
selectAtoms(bgf_file, selection, out_file):
    Extract some atoms from bgf_file. 
    IMPORTANT NOTE: Connection information are not returned.
            If you want to select atoms because of direct manipulation (i.e. adding 10 to xcoord)
            then you'd better use another function.
    """

    # variables
    myBGF = 0; myBGF2 = 0;

    # open bgf
    if isinstance(bgf_file, bgf.BgfFile):
        myBGF = bgf_file
    else:
        if not silent: print("Reading " + bgf_file + " ..")
        myBGF = bgf.BgfFile(bgf_file)

    # find 
    myBGF2 = bgf.BgfFile()
    executable = "if " + selection + ": " + "myBGF2.addAtom(atom)"
    for atom in myBGF.a:
        exec(executable)

    # if more than one atom is selected:
    # save
    if len(myBGF2) > 0:
        if isinstance(out_file, str):
            if not silent: print("Saving information to " + out_file + " ..")
            myBGF2.saveBGF(out_file)
            return 1;
        else:
            return myBGF2;
    else:
        print("No atoms are selected. Quit.")
        return 0;

    ### end of removeFragments


def charge(myBGF):
    charge = 0.0;
    for i in myBGF.a:
        charge += i.charge

    return float(charge)    # float


def mergeConnect(bgf1, ano1, bgf2, ano2):
    """
def mergeConnect(bgf1, ano1, bgf2, ano2):
    merge bgf1 and bgf2 & connect ano1 and ano2
    """

    atom1 = bgf1.getAtom(ano1)
    atom2 = bgf2.getAtom(ano2)

    # translate bgf2
    for atom in bgf2.a:
        atom.x += 2.0
        atom.y += 2.0
        atom.z += 2.0

    newbgf = bgf1.merge(bgf2)
    atom1.connect(atom2)

    return newbgf


def rotateBGF(bgf_file, atom1, atom2, vec, out_file, silent=False):
    """
rotateBGF: rotate v21 to the given vector (vec)
    """
 
    # open
    if isinstance(bgf_file, bgf.BgfFile):
        myBGF = bgf_file
    else:
        if not silent: print("reading " + bgf_file + " ..")
        myBGF = bgf.BgfFile(bgf_file)

    a1 = myBGF.a[myBGF.a2i[atom1]];
    a2 = myBGF.a[myBGF.a2i[atom2]];

    # original position of a1
    orig_x = a1.x; orig_y = a1.y; orig_z = a1.z

    v1 = (a1.x, a1.y, a1.z)

    ## move atom1 to origin
    #for atom in myBGF.a:
    #    atom.x -= v1[0]
    #    atom.y -= v1[1]
    #    atom.z -= v1[2]

    v21 = [ (a2.x - a1.x), (a2.y - a1.y), (a2.z - a1.z) ]

    # rotate v21 to x axis
    u1 = v21 / np.linalg.norm(v21)    ##

    v2 = [u1[1], -u1[0], 0]
    u2 = v2 / np.linalg.norm(v2)    ##

    v3 = np.cross(u1, u2)
    u3 = v3 / np.linalg.norm(v3)    ##

    U = np.array([[u1[0], u1[1], u1[2]], [u2[0], u2[1], u2[2]], [u3[0], u3[1], u3[2]]])

    # V = vec -> (1, 0, 0)
    # invV = inverse(V)
    V1 = vec / np.linalg.norm(vec)
    vec2 = [V1[1], -V1[0], 0]
    V2 = vec2 / np.linalg.norm(vec2)
    vec3 = np.cross(V1, V2)
    V3 = vec3 / np.linalg.norm(vec3)
    V = np.array([[V1[0], V1[1], V1[2]], [V2[0], V2[1], V2[2]], [V3[0], V3[1], V3[2]]])
    invV = np.linalg.inv(V)

    # rotate all atoms
    for atom in myBGF.a:
        a = np.matrix([atom.x, atom.y, atom.z]).T
        b = U*a
        c = invV*b
        atom.x = float(c[0])
        atom.y = float(c[1])
        atom.z = float(c[2])

    # move all atoms to the original position
    after_x = a1.x; after_y = a1.y; after_z = a1.z;
    delta_x = orig_x - after_x; delta_y = orig_y - after_y; delta_z = orig_z - after_z;
    for atom in myBGF.a:
        atom.x -= after_x
        atom.y -= after_y
        atom.z -= after_z

    ## save
    if isinstance(out_file, str):
        if not silent: print("saving information to " + out_file + " ..")
        myBGF.saveBGF(out_file)
        return 1;
    else:
        return myBGF;

    ### end of function


def make_periodic(*args):
    """
    MakePeriodic() returns 0
    MakePeriodic(bgf_file) returns BgfFile() object with CRYSTX {0 0 0 90 90 90}
    MakePeriodic(bgf_file, pbc) returns BgfFile() object with CRYSTX {x y z 90 90 90}
    MakePeriodic(bgf_file, CRYSTX) returns BgfFile() object with CRYSTX {x y z a b c}
    MakePeriodic(bgf_file, out_file, pbc) records out_file and returns 1
    MakePeriodic(bgf_file, out_file, CRYSTX) records out_file and returns 1
    """
    if len(args) == 0:
        nu.warn("No BGF file or instance specified.")
        return 0;
    elif len(args) >= 1:
        mybgf = args[0]
        if isinstance(mybgf, str):
            mybgf = bgf.BgfFile(mybgf)

        mybgf.PERIOD = "111"
        mybgf.AXES = "ZYX"
        mybgf.SGNAME = "P 1                  1    1"
        mybgf.CELLS = [-1, 1, -1, 1, -1, 1]

        if len(args) == 1:
            mybgf.CRYSTX = [0.0, 0.0, 0.0, 90.0, 90.0, 90.0]
            nu.warn("PBC set to [0.0, 0.0, 0.0, 90.0, 90.0, 90.0]")
            return mybgf
        elif len(args) >= 2:
            if len(args[-1]) == 3:
                mybgf.CRYSTX = args[1] + [90.0, 90.0, 90.0]
            elif len(args[-1]) == 6:
                mybgf.CRYSTX = args[1]
            else:
                nu.warn("Wrong pbc provided.")
                return 0;
            if len(args) == 2:
                return mybgf
            elif len(args) == 3:
                mybgf.saveBGF(args[1])
                return 1;
            else:
                nu.warn("Too many parameters (>4) provided.")
                return 0;


def copy_bgf_info(bgf_file, out_file=""):
    """
    Generates an empty BGF file with same information.
    copy_bgf_info(bgf_file): return new bgf file, DON'T SAVE
    copy_bgf_info(bgf_file, out_file="filename"): return new bgf file, DO SAVE
    """

    if isinstance(bgf_file, bgf.BgfFile):
        orig = bgf_file
    else:
        orig = bgf.BgfFile(bgf_file)
    new = bgf.BgfFile()

    new.BIOGRF = orig.BIOGRF
    new.DESCRP = orig.DESCRP
    new.REMARK = orig.REMARK
    new.FF = orig.FF
    new.FORMAT = orig.FORMAT
    new.PERIOD = orig.PERIOD
    new.AXES = orig.AXES
    new.SGNAME = orig.SGNAME
    new.CRYSTX = orig.CRYSTX
    new.OTHER = orig.OTHER
    new.CELLS = orig.CELLS

    if out_file:
        new.saveBGF(out_file)
        
    return new

            
def get_residue_names(bgf_file):
    """
    returns all residue names in the BGF file.
    """

    # bgf open
    if isinstance(bgf_file, bgf.BgfFile):
        mybgf = bgf_file
    else:
        mybgf = bgf.BgfFile(bgf_file)

    result = set()
    for atom in mybgf.a:
        result.add(atom.rName)

    return list(result)


def sort_atoms_in_residue_name(bgf_file, out_file=""):
    """
    Sort atoms in residue order in BGF file.
    Water atoms will be located in the last.

    sort_atoms_in_residue_name(bgf_file): 1) sort atoms 2) DO NOT SAVE 3) return BgfFile object
    sort_atoms_in_residue_name(bgf_file, out_file="filename"): 1) sort atoms 2) DO SAVE 3) return BgfFile object
    """
    
    if isinstance(bgf_file, bgf.BgfFile):
        orig = bgf_file
    else:
        orig = bgf.BgfFile(bgf_file)

    new = copy_bgf_info(bgf_file)

    print("Fetching residue names...")
    l_rNames = get_residue_names(bgf_file)
    print("\tfound residue names:" + str(l_rNames))
    l_rNames.sort()

    d_rNames_aNo = dict()   # d_rNames_aNo[rname] = [aNo]

    # make sure that residue name WAT is located on the last of the list
    print("Sorting residue names...")
    l_rNames_wat = ['WAT']
    for i in l_rNames_wat:
        try:
            l_rNames.remove(i)
        except ValueError:
            pass;
        else:
            l_rNames.append(i)
    print("\trefined residue names:" + str(l_rNames))

    # record aNos per rNames
    print("Reading atom information...")
    for r in l_rNames:
        temp = []
        for atom in orig.a:
            if r in atom.rName:
                temp.append(atom.aNo)
        
        d_rNames_aNo[r] = temp  # register

    # move atoms from orig to new
    print("Reordering atom information...")
    for r in l_rNames:
        l_ano = d_rNames_aNo[r]
        for i in l_ano:
            atom = orig.getAtom(i)
            new.addAtom(atom)

    new.renumber()

    # save if out_file provided
    if out_file:
        print("Saving BGF file to %s..." % out_file)
        new.saveBGF(out_file)
        return True

    return new


def remove_bad_contacts(bgf_file, out_file='', thresh=2.0, silent=False):

    import scipy.spatial.distance as dist

    # initialization
    if not isinstance(bgf_file, bgf.BgfFile):
        if not silent: print("Removing bad contacts from " + str(bgf_file) + " with distance threshold " + str(thresh) + " A and saving to " + str(out_file) + ".")
        myBGF = bgf.BgfFile(bgf_file);
    else:
        myBGF = bgf_file

    myBGF.renumber()
    myBGF = sort_atoms_in_residue_name(myBGF)  # water atoms should locate to the last

    # move coordinates to numpy array
    if not silent: print("Sorting atoms..")
    water = []; solute = []; all = []
    for atom in myBGF.a:
        all.append([atom.x, atom.y, atom.z])
        if 'OW' in atom.ffType or 'HW' in atom.ffType or 'WAT' in atom.rName:
            water.append([atom.x, atom.y, atom.z])
        else:
            solute.append([atom.x, atom.y, atom.z])

    # distance calculation
    if not silent: print("Calculating distances..")
    indices = np.where(dist.cdist(water, solute) <= thresh)[0]
    indices = np.unique(indices)    # water indices which is closer than thresh

    # paperworks to remove atoms from BGF
    anos = []
    del_list = []; del_rno_list = [];
    for i in indices:
        ano = i + 1 + len(solute)
        anos.append(ano)
    for i in anos:
        atom = myBGF.getAtom(i)
        del_rno_list += [atom.rNo]
        del_rno_list += atom.CONECT
    del_rno_list = np.unique(del_rno_list)
    del_rno_list = list(del_rno_list)
    
    if not del_rno_list:
        if not silent: print("There are no water molecules that corresponds to the criteria.")
        return 0;
    else:
        for atom in myBGF.a:
            if atom.rNo in del_rno_list:
                del_list.append(myBGF.a2i[atom.aNo])
        del_list.sort()
        del_list.reverse()

        if not silent: print("%d atoms will be removed from the BGF file." % len(np.unique(del_list)))
        myBGF.delAtoms(del_list, False)
        myBGF.renumber()

    myBGF.REMARK.append("remove_bad_contacts: threshold %5.1f" % thresh)
    # returns
    if not out_file:
        return myBGF;
    else:
        myBGF.saveBGF(out_file)
        if not silent: print("File saved to %s" % out_file)
        return 1;

    ### end of remove_bad_contacts


def stress_cell(bgf_file, ratio, out_file=''):
    '''
    scales the cell to the ratio.
    '''

    # initialization
    if not isinstance(bgf_file, bgf.BgfFile):
        myBGF = bgf.BgfFile(bgf_file);
    else:
        myBGF = bgf_file

    if not myBGF.CRYSTX:
        nu.die("No pbc information to stress the cell.")

    xs = ys = zs = 1.0
    parse = ratio.split()

    if len(parse) == 1:
        xs = ys = zs = float(parse[0])
    elif len(parse) == 3:
        xs = float(parse[0])
        ys = float(parse[1])
        zs = float(parse[2])
    else:
        nu.die("Error on specifying cell inflation ratio %s: must be 1 or 3" % ratio)

    for atom in myBGF.a:
        atom.x *= xs
        atom.y *= ys
        atom.z *= zs

    myBGF.CRYSTX[0] *= xs
    myBGF.CRYSTX[1] *= ys
    myBGF.CRYSTX[2] *= zs

    # returns
    if not out_file:
        return myBGF;
    else:
        myBGF.saveBGF(out_file)
        if not silent: print("File saved to %s" % out_file)
        return 1;


#def trjchunk2bgf(mybgf, chunk)



#-------------------------------------
#
# - Print information
#
if __name__ == '__main__':

    # get directory:
    directory = dir()

    # set imported stuff we don't want to see:
    imported = ['sys', 'bgf', 'bgf', 'string', 'dreiding', 'cu', 'os']

    # print __doc__ for the module:
    print("\n")
    print("-"*60)
    if 'version' not in directory:  version = '??????'
    print("%-45s%15s" % (os.path.basename(sys.argv[0]), 'ver: '+version))

    print("-"*60)
    print(__doc__)


    # import types:
    import types

    # create hash-table:
    hashtable = {}
    for item in directory:
        actual_item = eval(item)
        if item in imported:
            # don't show imported stuff:
            pass
        elif type(actual_item) is types.ModuleType:
            # don't discuss other modules:
            pass
        elif type(actual_item) is types.FunctionType:
            # show __doc__s for functions:
            hashtable[item] = actual_item.__doc__
        elif type(actual_item) is types.ClassType:
            # show __doc__s for classes:
            title = item+' class: '
            hashtable[item] = title +  ( '-' * (60-len(title)) )
            hashtable[item] += actual_item.__doc__

            # show __doc__s for class elements:
            for classItem in dir(actual_item):
                actual_class_item = eval(item+'.'+classItem)
                if type(actual_class_item) is types.ModuleType:
                    # don't discuss other modules:
                    pass
                elif type(actual_class_item) is types.UnboundMethodType \
                        or type(actual_class_item) is types.MethodType:
                    # show __doc__s for functions:
                    hashtable[item] += actual_class_item.__doc__ + '\n'

                elif classItem in ['__doc__','__module__']:
                    pass
                else:
                    # for other stuff show the value:
                    hashtable[item] += '\n'+classItem+' = '+str(actual_class_item)+'\n'

            hashtable[item] +=  ( '-'*60 )+'\n\n'
            
        elif item[0] != '_':
            # for other stuff show the value:
            hashtable[item] = '\n'+item+' = '+str(actual_item)+'\n'

    # print info out
    keys = hashtable.keys()
    keys.sort()
    
    print("Contents:")
    print("-"*60 )
    for item in keys:
        print(hashtable[item])

    print("\n")
    print("-"*60 )
    print("contact: noische@kaist.ac.kr")
    print("-"*60 )
    print("\n")
    
    # done!

