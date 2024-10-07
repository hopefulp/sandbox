#### to modify vasp input file
import sys
import vasp_job
import re
import os
from common import whereami
'''
def get_ntatom(st):
def parse_poscar(pos, opt):
def get_iatoms_in_group(zmin, zmax, coord, loc):
def obtain_atomlist0(zminmax, poscar, atom_species, loc):
def get_atoms_poscar(line):
def get_poscar(poscar, job='new', sub=0):
def get_poscar0(poscar):
def pos2dirname(poscar):
def get_dnames4pos(poscars):
def fixedMD_POSCAR(poscar, atom, atoms=None):
'''
def get_ntatom(st):
    natom_list = list(map(int, st.strip().split()))
    ntotal = sum(natom_list)
    return ntotal, natom_list

### read_poscar was changed into parse_poscar
def parse_poscar(pos, opt):
    '''
    extract lattice vectors, pre-part, atom_list, coordinates
    opt    pre, coord,
    '''
    with open(pos, "r") as f:
        lines = f.readlines()
        p_axes=[]   # principal axis

        for i in range(2,5):
            paxis = list(map(float,lines[i].strip().split()))
            p_axes.append(paxis)

        if opt == 'paxes':
            return p_axes
        else:
            for i in [5, 0]:
                if not any(s.isdigit() for s in lines[i]):
                    atom_list = lines[i].strip().split()
                    break
            for i in [5, 6]:
                if any(s.isdigit() for s in lines[i]):
                    ntatom, natom_list = get_ntatom(lines[i])
                    break
            for i in [8, 9]:
                if any(s.isdigit() for s in lines[i]):
                    iend_pre = i-1
                    break
            ### if direct, use z-size of paxis_c to scale down the z_coord of cartesian from visual program: ase gui
            if re.match('D', lines[iend_pre], re.I):
                p_z = p_axes[2][2]
            else:
                p_z = 1.0
            pre = lines[:iend_pre+1]
            coord = lines[iend_pre+1:iend_pre+ntatom+1]

            if opt == 'pre':
                return pre
            elif opt == 'coord':
                return coord
            elif opt == 'atomcoord':
                return atom_list, natom_list, coord, p_z
            elif opt == 'alist':
                return atom_list, natom_list

def get_iatoms_in_group(zmin, zmax, coord, loc):
    ind=[]
    for i, xyzs in enumerate(coord):
        line_ele = xyzs.strip().split()
        #print(f"{line_ele}")
        if loc == 'in':
            if zmin < float(line_ele[2]) and float(line_ele[2]) < zmax:
                ind.append(i)
        elif loc == 'out':
            if float(line_ele[2]) < zmin or zmax < float(line_ele[2]):
                ind.append(i)
    return ind


def obtain_atomlist0(zminmax, poscar, atom_species, loc):
    '''
    obtain_atomlist0: '0' denotes starting atom index
    input
        read poscar
    return
        atom list inbetween zmin & zmax: index from 0 ~
        principal axes    
    '''
    if len(zminmax) == 2:
        zmin, zmax = (zminmax[0], zminmax[1])
    else:
        print("z-coord error: {zminmax}, input two z-values with -z ")
        return 1
    atom_list, natom_list, coord, direct_z = parse_poscar('atomcoord', poscar)
    ### in case direct coordinates in POSCAR, reduce cartisian by /direct_z
    if direct_z != 1.0:
        zmin /= direct_z
        zmax /= direct_z
        print(f"{whereami():>15}(): zmin, max = {zmin} {zmax} in direct coordinates in POSCAR: in {__name__}.py")
    iatom=0
    ind0_select=[]
    print(f"{whereami():>15}(): POSCAR {atom_list} {natom_list} {len(coord)} coordinates: in {__name__}.py ")
    for i, atom in enumerate(atom_list):
        if atom in atom_species:
            ind0_group = get_iatoms_in_group(zmin, zmax, coord[iatom:iatom+natom_list[i]], loc)
            ind0_select.extend([x+iatom for x in ind0_group])
        iatom += natom_list[i]

    print(f"{whereami():>15}(): indices {ind0_select} total {len(ind0_select)}")
    return ind0_select

     

def get_atoms_poscar(line):
    if atom in lines[0]:
        alist = lines[0].split()
    return alist


def get_poscar(poscar, job='new', sub=0):
    '''
    copy {poscar} to POSCAR at cwd
    poscar  input of POSCAR.name
        if dir:
            read dir/POSCAR or read dir/CONTCAR
    return new dirname
    '''
    dname=None
    ### 
    if not poscar :
        print("POSCAR will be used")
    ### if poscar is dir
    else:
        #print(f"{__name__}: {poscar}")
        if os.path.isdir(poscar):
            if sub == 0:
                nposcar = f'{poscar}/POSCAR'
            else:
                nposcar = f'{poscar}/CONTCAR'
            dname = poscar + job
        else:
            nposcar = poscar
            dname = pos2dirname(poscar)
        comm = 'cp %s POSCAR' % nposcar
        print(f"{__name__}:{whereami()}:: {comm}")
        os.system(comm)
        print(f'{__name__}:{whereami()}:: POSCAR was made from {nposcar}')
    # confirm POSCAR is made        
    if not os.access('POSCAR', os.F_OK):
        print('POSCAR is not here')
        exit(21)
    return dname

def get_poscar0(poscar):
    '''
    copy {poscar} to POSCAR at cwd
    '''
    ### confirm file location
    ### if poscar is Null, use POSCAR
    if not poscar :
        print("POSCAR will be used")
    elif not os.access('%s' % poscar, os.F_OK):
        print('poscar is not detectable')
        exit(2)
    else:
        comm = 'cp %s POSCAR' % poscar
        os.system(comm)
        print('POSCAR is made')
    # confirm POSCAR is made        
    if not os.access('POSCAR', os.F_OK):
        print('POSCAR is not here')
        exit(21)

    return 0

def pos2dirname(poscar):
    ### obtain dirname from POSCAR.dir
    if re.match("POSCAR", poscar):
        if len(poscar) == 6:
            dirname='pos'
        else:
            dirname = poscar[7:]
    else:
        dirname = poscar
    return dirname

def get_dnames4pos(poscars):
    dnames=[]
    for poscar in poscars:
        dnames.append(pos2dirname(poscar))
    return dnames


def fixedMD_POSCAR(poscar, atom, atoms=None):
def fixedMD_POSCAR(poscar, atom, atoms=None):

    '''
    poscar  to be modified
    atom    kind to be moved for zpe calculation
    '''
    with open(poscar, 'r') as f:
        lines = f.readlines()
    with open('POSCAR', 'w') as f:
        iatom = 0
        for i, line in enumerate(lines):
            if i == 0:
                f.write(line)
                continue
            elif i < 5:
                f.write(line)
            ### atom lines could exist or not
            elif i==5 and any(j.isalpha() for j in line):
                atoms = line.strip().split()
                f.write(line)
            elif i <= 6 and any(j.isdigit() for j in line):
                natoms = list(map(int, line.strip().split()))
                if len(atoms) != len(natoms):
                    sys.exit(0)
                f.write(line)
                ### calculate
                ntotal = sum(natoms)
                ind = atoms.index(atom)
                nzpe = natoms[ind]
                isum = 0
                for i, na in enumerate(natoms):
                    if i < ind:
                        isum += na 
                f.write("Selective dynamics\n")
                print(f"ind {ind} pre sum {isum} ")
            elif i <= 7 and any(i.isalpha() for i in line):
                f.write(line)
            else:
                if iatom < isum:
                    new_line = line.rstrip() + " F F F\n"
                    f.write(new_line)
                elif iatom < isum + nzpe:
                    new_line = line.rstrip() + " T T T\n"
                    f.write(new_line)
                elif iatom < ntotal:
                    new_line = line.rstrip() + " F F F\n"
                    f.write(new_line)
                else:
                    break
                iatom += 1
    return 0
