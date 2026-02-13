import numpy as np
from common import whereami
import sys

### DOSCAR ordering

orbital = ['s', 'p', 'd', 'f']
ltos = {'0':'s', '1':'p', '2':'d', '3':'f'}
l_map = {'s': [1], 'p': [2,3,4], 'd': [5,6,7,8,9]}

#mag_orbital = [['s'],['y','z','x'],['xy', 'yz', 'z2', 'xz', 'x2-y2']] 
mag_orbital = ['s','y','z','x','xy', 'yz', 'z2', 'xz', 'x2-y2'] 

is_doscar = [0]
ip_doscar = [1,2,3]
id_doscar = [4,5,6,7,8]

n_header = 5
dos_err = 100

def doscar_Bheader(lines8):
    '''
    read DOSCAR header + 2 line
        5 headline
        6 Block headline
        7-8 1-2 TDOS values to check error

    output  int:    natom, ngrid
            float:  Emax, Emin, Ef
            str:    bheadline (block headline)
            boolean: check doserr exists
    '''
    #print(f"header in {whereami()}(): {lines8}")
    if not lines8:
        print(f"proper DOSCAR is not provided, move to dos directory")
        sys.exit(111)
    for i, line in enumerate(lines8):
        if i == 0:
            natom = int(line.strip().split()[0])
        ### DOS Block headline
        elif i == n_header: 
            info = line.strip().split()
            Emax    = float(info[0])
            Emin    = float(info[1])
            ngrid   = int(info[2])
            Ef      = float(info[3])
        ### check whether DOS at first energy has error
        elif i == n_header +1:
            elist   = line.strip().split()
            dos1    = float(elist[1])
            ncol_tdos = len(elist)
        elif i == n_header + 2:
            dos2    = float(line.strip().split()[1])
    if dos2 * dos_err < dos1:
        Ldos_err = True
    else:
        Ldos_err = False
    if ncol_tdos == 3:
        Lspin = False
    elif ncol_tdos == 5:
        Lspin = True
    else:
        print(f"DOSCAR parsing error: ncol of tdos {ncol_tdos}")
        sys.exit(100)
    #print(f"{whereami():>15}(): DOS err {Ldos_err} at 1st energy {dos1} then {dos2} in {__name__}.py")
    return natom, Emax, Emin, ngrid, Ef, Ldos_err, Lspin

def change_Bheadline(old, Emin, Erep, ngrid, new_ngrid):
    print(f"{old}")
    new_grid = old.replace(ngrid, new_ngrid)
    print(f"{new_grid}")
    new_headline = new_grid.replace(Emin, Erep)
    print(f"{new_headline} with {Erep} for {Emin} ")
    return new_headline  


def read_doscar(filename="DOSCAR", atom_indices=None, l=None, option='dos'):
    """
    option: head    
            TDOS
            PLDOS
    Read and sum projected DOS for selected atoms and angular momentum.
    
    Parameters:
        filename: str – path to DOSCAR
        atom_indices: list[int] – list of atom indices (1-based, like in VASP)
        l: str – orbital type ('s', 'p', or 'd')
        Lspin   .T. or .F.
    
    Returns:
        dos_data: np.ndarray – array with shape (n_ene, 2), columns = [energy, summed_dos]
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    ### Part 1: read header
    ### head information is needed without l
    # Read header: 1st and n_header line; 3 more lines for checking DOS value error
    natom, Emax, Emin, n_ene, Ef, Ldoserr, Lspin = doscar_Bheader(lines[0:n_header+3])

    #if max(atom_indices) > natom:
    #    print(f"ERROR: max atom index {max(atom_indices)} exceeds natom {natom}")
    #    sys.exit(1)

    if option == 'head':
        return Ef, Lspin, Ldoserr

    # Read energy grid from total DOS block
    tdos_start = n_header
    tdos_end = tdos_start + 1 + n_ene
    tdos_block = lines[tdos_start + 1: tdos_end]
    energy = np.loadtxt(tdos_block, usecols=0)
    summed_dos = np.zeros(n_ene)

    if atom_indices[0] == 0:
        if not Lspin:
            summed_dos = np.loadtxt(tdos_block, usecols=[1])
        else:
            t_dos = np.loadtxt(tdos_block, usecols=[1,3])
            summed_dos = t_dos.sum(axis=1)

        dos_data = np.column_stack((energy, summed_dos))
        print(f"TDOS: shape {dos_data.shape}, first row: {dos_data[0]}")
        return dos_data
    
    if l.isdigit():
        l = ltos[l]

    #l_map = {'s': [1], 'p': [2,3,4], 'd': [5,6,7,8,9]}
    ### DOSCAR read 9 columns but PROCAR has 10-th for total
    if l_map.get(l):
        col_idx = l_map[l]
    else:
        if l == 't':
            l = 'spd'
        col_idx = []
        if 's' in l:
            col_idx.extend(l_map['s'])
        if 'p' in l:
            col_idx.extend(l_map['p'])
        if 'd' in l:
            col_idx.extend(l_map['d'])
    
    print(f"column indices {col_idx}")

    #print(f"energy {energy[0]}, {energy[-1]}")
    # Prepare array for summed DOS
    #summed_dos = np.zeros(n_ene)
    
    # Start index for PDOS blocks
    pdos_start = n_header + 1 + n_ene   # first line block head

    for i in atom_indices:
        # VASP uses 1-based indexing; convert to 0-based
        block_start = pdos_start + (i - 1) * (n_ene + 1)
        atom_block = lines[block_start + 1: block_start + 1 + n_ene]

        pdos = np.loadtxt(atom_block, usecols=col_idx)
        if pdos.ndim == 1:
            summed_dos += pdos
        else:
            summed_dos += pdos.sum(axis=1)
    print(f"shape of energy {energy.shape} summed_dos {summed_dos.shape} in {whereami()}() in module {__name__}")
    dos_data = np.column_stack((energy, summed_dos))
    print(f"dos_data {dos_data.shape} in {whereami()}() in module {__name__}")
    return dos_data

### from mod_doscar.py
def obtain_doscar_head(fname):
    '''
    read DOSCAR
    return  int:    natom, ngrid
            float:  Emax, Emin, Ef
            str:    bheadline (block headline)
            boolean: check doserr exists
    '''
    with open(fname, 'r') as f:
        ### analyze DOSCAR
        for i, line in enumerate(f):
            if i < nheadline:
                ### natom natom 1 0
                if i == 0:
                    natom = int(line.strip().split()[0])
                ### : if not just fout.write(line)
            elif i == nheadline:
                bheadline = line
                info = line.strip().split()
                Emax    = float(info[0])
                Emin    = float(info[1])
                ngrid   = int(info[2])
                Ef      = float(info[3])
            ### check whether DOS at first energy has error
            elif i == nheadline +1:
                elist   = line.strip().split()
                dos1    = float(elist[1])
                ncol_tdos = len(elist)
            elif i == nheadline + 2:
                dos2    = float(line.strip().split()[1])
                break
    if dos2 * dos_err < dos1:
        Ldos_err = True
    else:
        Ldos_err = False
    if ncol_tdos == 3:
        Lspin = False
    elif ncol_tdos == 5:
        Lspin = True
    else:
        print(f"DOSCAR parsing error: ncol of tdos {ncol_tdos}")
        sys.exit(100)
    print(f"{whereami():>15}(): DOS err {Ldos_err} at 1st energy {dos1} then {dos2} in {__name__}.py")
    return natom, Emax, Emin, ngrid, Ef, bheadline, Ldos_err, Lspin

def change_Bheadline(old, Emin, Erep, ngrid, new_ngrid):
    print(f"{old}")
    new_grid = old.replace(ngrid, new_ngrid)
    print(f"{new_grid}")
    new_headline = new_grid.replace(Emin, Erep)
    print(f"{new_headline} with {Erep} for {Emin} ")
    return new_headline  
