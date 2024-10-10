#### to modify vasp input file
import sys, re, os
import vasp_job
import numpy as np
from common import whereami
import fortranformat as ff
from univ_const import k_B, amu     # k_B (J/K), amu = atomic mass unit (kg) will be multiplied by atomic mass
from ase.data import atomic_masses, chemical_symbols    # atomic weight, symbol with same order in the lists
from my_stat import get_MBD_1D
from parsing import startnum

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

def get_zmax(poscar):
    coords = parse_poscar(poscar, block='coord'):
    
    


### read_poscar was changed into parse_poscar
def parse_poscar(pos, block = None, opt=None):
    '''
    old
        extract lattice vectors, pre-part, atom_list, coordinates
        opt    pre, coord,
    new
        parse by block
        returns line of string
    opt returns values
    '''
    with open(pos, "r") as f:
        lines = f.readlines()

        ### whether atom list line exists
        Latomline = 0                   # normally exists -> do not change line index
        if not all(x.isalpha() or x.isspace() for x in lines[5].strip()):
            Latomline = -1              # if not, lower index by 1

        ### ntotal is required to extract coord block
        i = 6 + Latomline
        natoms = list(map(int, lines[i].strip().split()))
        ntotal = sum(natoms)

        ### whether selective MD line exists
        Lselect = 0                 # normally doesnot exist
        for i in [6, 7]:
            if re.match('s', lines[i], re.I):
                Lselect = 1         # if yes, increase index by 1
        
        p_axes=[]   # principal axis
        if block == 'title':
            atoms = lines[0].strip().split()
            return lines[0], atoms
        elif block == 'scale':
            scale = lines[1].strip()
            return lines[1], scale
        elif block == 'paxes':
            nline = []
            for i in range(2,5):
                nline.append(lines[i])
                paxis = list(map(float,lines[i].strip().split()))
                p_axes.append(paxis)
            return nline, p_axes 
        elif block == 'atoms':
            if Latomline == 0:
                atoms = lines[5].strip().split()
                return lines[5], atoms
            else:
                return None, None
        elif block == 'natoms':
            iline = 6 + Latomline
            natoms = list(map(int, lines[iline].strip().split()))
            return lines[iline], natoms
        elif block == 'cd':
            iline = 7 + Latomline + Lselect
            if not re.match('[D|C]', lines[iline], re.I):       # check for logic
                print(f'line counting error in {whereami()}() in {__name__}.py')
                sys.exit(100)   
            cd = lines[iline].strip()
            return lines[iline], cd
        elif block == 'coord':
            istart = 8 + Latomline + Lselect
            coord = lines[istart:istart+ntotal]     # this is string of line with \n
            ### return 2D(natom by 3) coordinates
            if opt:
                print(f'change into numbers in {whereami()}()')
                coordnum = []
                for line in coord:
                    lxyz = line.strip().split()
                    coordnum.append(lxyz)
                coord=coordnum
            return coord, ntotal                    # return lines

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
    job     mode       s/a  select or add (job has 'add')
                zpe     s   fixedMD or selective MD for selective dynamics for ZPE calculation
                bomb    s   velocity to bombardment experiment to -z axis
                addbomb a   add molecule and velocity to -z axis
    matoms       aAtomN for add atom, sAtomN for select atom in POSCAR
                a   find lattice constand and distribute added atoms
                    append N atoms
                s   select from atoms, natoms list in POSCAR
                    Hf, O, Mo, S, O -> O1, O2 
                    atomlist in POSCAR for movable in zpe calculation
                velocity vel depending on T
    option      job = bomb addbomb:  T
    iatom       atom index
    atoms       atom list in POSCAR
    natoms      number of atoms list in POSCAR
    '''
    print(f"Write to {outf} in {whereami()}()")

    lines = []
    #with open(outf, 'w') as f :
    ### 1st line
    ### obtain each block and write from parse_POSCAR
    line1, atoms1 = parse_poscar(poscar, block='title'); lines.append(line1) # atoms might appear 1st line
    line, scale = parse_poscar(poscar, block='scale'); lines.append(line)
    nline, paxes = parse_poscar(poscar, block='paxes'); lines.extend(nline)
    line, atoms = parse_poscar(poscar, block='atoms')
    if not line:
        line = line1
        atoms = atoms1      # line1 is saved for atom line
    ### add atom
    if 'add' in job:
        ### split alphabet and numeric: O7, O29, He7
        i = startnum(matoms)
        addatom = matoms[:i]
        naddatom = matoms[i:]
        line += f' {addatom}' + '\n'
        atomsold = atoms[:]
        atoms.append(addatom)
    lines.append(line)
    line, natoms = parse_poscar(poscar, block='natoms')
    ntotalold = sum(natoms)
    if 'add' in job:
        line += f' {naddatom}' + '\n'
        natomsold = natoms[:]
        natoms.append(naddatom)
    lines.append(line)
    ### for sigma for MBD (Maxwell-Boltzmann distribution)
    atom_oldfull_list = make_atomfullist(atomsold, natomsold)       # As for original POSCAR,

    ### for selection mode
    if not 'add' in job:
        ### calculate index for selected atoms: movable in zpe, velocity in bombardment
        ind, nselatoms = select_atoms_poscar(atoms, natoms, matoms)
        nselatoms = natoms[ind]              # natoms to be moved for zpe
        npre_unsel = 0
        for i, na in enumerate(natoms):
            if i < ind:
                npre_unsel += na                  # npre_unsel = natoms before selected (movable) atoms
        if job == 'zpe':
            lines.append("Selective dynamics\n")
        print(f"ind {ind}  {npre_unsel} in {whereami()}()")
    
    ### for Cartesian or Direct    
    line, _ = parse_poscar(poscar, block='cd'); lines.append(line)

    ### looping in atomic position block
    coords, _ = parse_poscar(poscar, block='coord'); 
    
    ## modify selective dynamics for ZPE
    if job == 'zpe':
        new_coords=[]
        for i, line_coord in enumerate(coords):
            if i < npre_unsel:
                new_line = line_coord.rstrip() + " F F F\n"
            elif i < npre_unsel + nselatoms:
                new_line = line_coord.rstrip() + " T T T\n"
            elif i < ntotalold:
                new_line = line_coord.rstrip() + " F F F\n"
            new_coords.append(new_line)
        coords = new_coords
    lines.extend(coords)

    ### Make additional coordinates for added atoms
    if 'add' in job:
        ### for Orthorhombic crystal structure
        zcoord = 5.0    # bombing atoms to z-axis from surface
        z_surf = obtain_zmax(poscar)
        add_coords = surf_distributtion(atomname, n, paxes)

        lines.extend(add_coords)
        
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
    poscar selective bombarding is done -> upgrade to modify_POSCAR
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
