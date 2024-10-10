#### to modify vasp input file
import sys, re, os
import vasp_job
import numpy as np
from common import whereami
import fortranformat as ff
from univ_const import k_B, amu     # k_B (J/K), amu = atomic mass unit (kg) will be multiplied by atomic mass
from ase.data import atomic_masses, chemical_symbols    # atomic weight, symbol with same order in the lists
from my_stat import get_MBD_1D
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
def modify_POSCAR <- fixedMD_POSCAR(poscar, atom, atoms=None):
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
        if opt == 'title':
            atoms = lines[0].strip().split()
            return lines[0], atoms
        
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


def make_atomfullist(atom_list, natom_list):
    '''
    from atomlist and natomlist in POSCAR
    make full atom list
    '''
    atom_fullist=[]
    if len(atom_list) == len(natom_list):
        for natom, symbol in zip(natom_list, atom_list):
            for j in range(natom):
                atom_fullist.append(symbol)
        return atom_fullist
    else:
        print(f"Parsing error:: different nspecies {len(atom_list)} != list of natoms {len(natom_list)}")
        sys.exit(100)

            

def select_atoms_poscar(atom_list, natom_list, sel_atom):
    '''
    from POSCAR atom and natom list return atom index and natoms for selection
    atom_list   in POSCAR
    natom_list  in POSCAR
    sel_atom    index or atom name with counting number
    '''
    if sel_atom.isalpha():          # O, Hf, etc
        ind         = atom_list.index(sel_atom) # atomlist is atom kind to be moved (T T T for ZPE)
    elif sel_atom.isdigit():        # 0, 1, ...
        ind         = sel_atom
    else:                           # O1, O2, ...
        atom        = sel_atom[:-1]
        order       = sel_atom[-1]
        dicatoms    = { n: rep[n] for rep in [{}] for i, n in enumerate(atom_list) if rep.setdefault(n, []).append(i) or len(rep[n])==2}
        print(dicatoms)
        ind         = dicatoms[atom][order]
    nselatoms   = natom_list[ind]          # natoms to be moved for zpe
    print(f"{ind} : {nselatoms} selected in {whereami()}()")
    return ind, nselatoms
        


#def fixedMD_POSCAR(poscar, atom, atoms=None?):
def modify_POSCAR(poscar, job='zpe', matoms=None, outf='POSCAR', option=None):
    '''
    Modularize POSCAR part
    poscar      to be modified
    job         zpe     fixedMD or selective MD for selective dynamics for ZPE calculation
                bomb    velocity to bombardment experiment to -z axis
                addbomb add molecule and velocity to -z axis
    matoms       aAtomN for add atom, sAtomN for select atom in POSCAR
                a   find lattice constand and distribute added atoms
                    append N atoms
                s   select from atoms, natoms list in POSCAR
                    Hf, O, Mo, S, O -> O1, O2 
                    atomlist in POSCAR for movable in zpe calculation
                velocity vel depending on T
    option      job = vel, addbomb  T
    iatom       atom index
    atoms       atom list in POSCAR
    natoms      number of atoms list in POSCAR
    '''
    print(f"Write to {outf} in {whereami()}()")

    if re.match('a') in atoms:
        add_atoms = matoms[1:]

        ajob = 'add'
    elif re.match('s') in atoms:
        sel_atoms = matoms[1:]
        ajob = 'select'

    with open(outf, 'w') as f :
        ### obtain each block and write from parse_POSCAR
        ### 1st line
        line, atoms = parse_poscar(poscar, block='title')
            if i == 0:
                f.write(line)
                ### this might be atom list in POSCAR
                atoms = line.strip().split()
                continue
            elif i < 5:
                f.write(line)
            ### atom lines could exist or not
            elif i==5 and all(x.isalpha() or x.isspace() for x in line.strip()):
                ### if not here, 1st line is used for atoms
                atoms = line.strip().split()
                if 'add' in job:
                    line = line.rstrip() + f' {add_atoms[:-1]}' + '\n' 
                f.write(line)
            elif i <= 6 and all(x.isdigit() or x.isspace() for x in line.strip()):
                natoms = list(map(int, line.strip().split()))
                if len(atoms) != len(natoms):
                    print(f"Error:: len(atoms) {len(atoms)} != {len(natoms)} len(natoms)")
                    sys.exit(100)
                #if 'add' in job:
                #    line = line.rstrip() + f' {
                f.write(line)
                ### calculate index for selected atoms: movable in zpe, velocity in bombardment
                atom_flist = make_atomfullist(atoms, natoms)  # for sigma for MBD (Maxwell-Boltzmann distribution)
                ntotal = sum(natoms)
                ind, nselatoms = select_atoms_poscar(atoms, natoms, sel_atom)
                nselatoms = natoms[ind]              # natoms to be moved for zpe
                npre_unsel = 0
                for i, na in enumerate(natoms):
                    if i < ind:
                        npre_unsel += na                  # npre_unsel = natoms before selected (movable) atoms
                if job == 'zpe':
                    f.write("Selective dynamics\n")
                print(f"ind {ind}  {npre_unsel} in {whereami()}()")
            ### for Cartesian or Direct    
            elif i <= 7 and any(i.isalpha() for i in line):
                f.write(line)
            ### looping in atomic position block
            else:
                if job == 'zpe':
                    if iatom < npre_unsel:
                        new_line = line.rstrip() + " F F F\n"
                    elif iatom < npre_unsel + nselatoms:
                        new_line = line.rstrip() + " T T T\n"
                    elif iatom < ntotal:
                        new_line = line.rstrip() + " F F F\n"
                    else:
                        break
                    f.write(new_line)
                if 'add' in job:
                    paxes2D = parse_poscar(poscar, paxes)

                if 'bomb' in job:
                    f.write(line)
                iatom += 1
                if iatom == ntotal and 'vel' in job:
                    ### might need addatoms
                    ### addtional velocity block as cartisian coordinate (A/fs)
                    f.write("\n")
                    T = option
                    mu = 0.0
                    amukT = amu/(k_B * T)
                    ms2angfs = 1.E-5         # from m/s to (1E10/1E15) Ang/fs
                    lineformat = ff.FortranRecordWriter('3E16.8')
                    for i, atom in enumerate(atom_flist):
                        ### mass is required: depending on mass, sigma=1/(m/k_B*T) is changed
                        atomic_weight = atomic_masses[chemical_symbols.index(atom)]
                        sigma = 1./np.sqrt(atomic_weight*amukT) * ms2angfs 
                        vx, vy, vz = get_MBD_1D(loc=mu, scale=sigma, size=1)    # N.B. each v's are list of size 
                        ### only selected atoms have vel = -z direction for bombardment
                        if i < npre_unsel:
                            s = lineformat.write([vx[0], vy[0], vz[0]]) + "\n"
                        elif i < npre_unsel + nselatoms:
                            v = np.sqrt(vx[0]**2 + vy[0]**2 + vz[0]**2)
                            #s = lineformat.write([0.0, 0.0, -vz[0])]) + "\n"
                            s = lineformat.write([0.0, 0.0, -v]) + "\n"
                        elif iatom < ntotal:
                            s = lineformat.write([vx[0], vy[0], vz[0]]) + "\n"
                        else:
                            break
                        f.write(s)
                    break
    return 0



#def fixedMD_POSCAR(poscar, atom, atoms=None?):
def modify_POSCAR2(poscar, job='zpe', atoms=None, outf='POSCAR', option=None):
    '''
    poscar      to be modified
    job         zpe     fixedMD or selective MD for selective dynamics for ZPE calculation
                bomb    velocity to bombardment experiment to -z axis
                addbomb add molecule and velocity to -z axis
    atoms       aAtomN for add atom, sAtomN for select atom in POSCAR
                a   find lattice constand and distribute added atoms
                    append N atoms
                s   select from atoms, natoms list in POSCAR
                    Hf, O, Mo, S, O -> O1, O2 
                    atomlist in POSCAR for movable in zpe calculation
                velocity vel depending on T
    option      job = vel, addbomb  T
    iatom       atom index
    atoms       atom list in POSCAR
    natoms      number of atoms list in POSCAR
    '''
    print(f"Write to {outf} in {whereami()}()")

    if re.match('a') in atoms:
        add_atoms = atoms[1:]

        ajob = 'add'
    elif re.match('s') in atoms:
        sel_atoms = atoms[1:]
        ajob = 'select'

    with open(poscar, 'r') as rf, open(outf, 'w') as f :
        lines = rf.readlines()
        iatom = 0
        for i, line in enumerate(lines):
            ### read poscar line by line
            if i == 0:
                f.write(line)
                ### this might be atom list in POSCAR
                atoms = line.strip().split()
                continue
            elif i < 5:
                f.write(line)
            ### atom lines could exist or not
            elif i==5 and all(x.isalpha() or x.isspace() for x in line.strip()):
                ### if not here, 1st line is used for atoms
                atoms = line.strip().split()
                if 'add' in job:
                    line = line.rstrip() + f' {add_atoms[:-1]}' + '\n' 
                f.write(line)
            elif i <= 6 and all(x.isdigit() or x.isspace() for x in line.strip()):
                natoms = list(map(int, line.strip().split()))
                if len(atoms) != len(natoms):
                    print(f"Error:: len(atoms) {len(atoms)} != {len(natoms)} len(natoms)")
                    sys.exit(100)
                #if 'add' in job:
                #    line = line.rstrip() + f' {
                f.write(line)
                ### calculate index for selected atoms: movable in zpe, velocity in bombardment
                atom_flist = make_atomfullist(atoms, natoms)  # for sigma for MBD (Maxwell-Boltzmann distribution)
                ntotal = sum(natoms)
                ind, nselatoms = select_atoms_poscar(atoms, natoms, sel_atom)
                nselatoms = natoms[ind]              # natoms to be moved for zpe
                npre_unsel = 0
                for i, na in enumerate(natoms):
                    if i < ind:
                        npre_unsel += na                  # npre_unsel = natoms before selected (movable) atoms
                if job == 'zpe':
                    f.write("Selective dynamics\n")
                print(f"ind {ind}  {npre_unsel} in {whereami()}()")
            ### for Cartesian or Direct    
            elif i <= 7 and any(i.isalpha() for i in line):
                f.write(line)
            ### looping in atomic position block
            else:
                if job == 'zpe':
                    if iatom < npre_unsel:
                        new_line = line.rstrip() + " F F F\n"
                    elif iatom < npre_unsel + nselatoms:
                        new_line = line.rstrip() + " T T T\n"
                    elif iatom < ntotal:
                        new_line = line.rstrip() + " F F F\n"
                    else:
                        break
                    f.write(new_line)
                if 'add' in job:
                    paxes2D = parse_poscar(poscar, paxes)

                if 'bomb' in job:
                    f.write(line)
                iatom += 1
                if iatom == ntotal and 'vel' in job:
                    ### might need addatoms
                    ### addtional velocity block as cartisian coordinate (A/fs)
                    f.write("\n")
                    T = option
                    mu = 0.0
                    amukT = amu/(k_B * T)
                    ms2angfs = 1.E-5         # from m/s to (1E10/1E15) Ang/fs
                    lineformat = ff.FortranRecordWriter('3E16.8')
                    for i, atom in enumerate(atom_flist):
                        ### mass is required: depending on mass, sigma=1/(m/k_B*T) is changed
                        atomic_weight = atomic_masses[chemical_symbols.index(atom)]
                        sigma = 1./np.sqrt(atomic_weight*amukT) * ms2angfs 
                        vx, vy, vz = get_MBD_1D(loc=mu, scale=sigma, size=1)    # N.B. each v's are list of size 
                        ### only selected atoms have vel = -z direction for bombardment
                        if i < npre_unsel:
                            s = lineformat.write([vx[0], vy[0], vz[0]]) + "\n"
                        elif i < npre_unsel + nselatoms:
                            v = np.sqrt(vx[0]**2 + vy[0]**2 + vz[0]**2)
                            #s = lineformat.write([0.0, 0.0, -vz[0])]) + "\n"
                            s = lineformat.write([0.0, 0.0, -v]) + "\n"
                        elif iatom < ntotal:
                            s = lineformat.write([vx[0], vy[0], vz[0]]) + "\n"
                        else:
                            break
                        f.write(s)
                    break
    return 0
