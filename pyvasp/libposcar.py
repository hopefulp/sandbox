#### to modify vasp input file
import sys, re, os
import vasp_job
import numpy as np
from common import whereami
import fortranformat as ff
from univ_const import k_B, amu     # k_B (J/K), amu = atomic mass unit (kg) will be multiplied by atomic mass
from ase.data import atomic_masses, chemical_symbols, vdw_radii    # atomic weight, symbol with same order in the lists
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



### read_poscar was changed into parse_poscar
def parse_poscar(pos, block = None, opt=None):
    '''
    old
        extract lattice vectors, pre-part, atom_list, coordinates
        opt    pre, coord,
    new
        parse by block
        block: pre: title - natoms
            title
            scale
            paxes
            atoms
            natoms
            cd      returns 'C' or 'D'
            coord   natoms line for coord_x, _y, _z
            vel     natoms line for vel (A/fs) not equipped yet
    '''
    with open(pos, "r") as f:
        lines = f.readlines()

        ### whether atom list line exists
        Latomline = 0                   # in std, exists -> do not change line index
        if not all(x.isalpha() or x.isspace() for x in lines[5].strip()):
            Latomline = -1              # if not, lower index by 1

        ### ntotal in natom line is required to extract coord block
        inatom = 6 + Latomline
        natoms = list(map(int, lines[inatom].strip().split()))
        ntotal = sum(natoms)

        ### for pre-part return here
        if block == 'pre':
            return lines[:inatom+1]

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
            cd = lines[iline].strip()[0].capitalize()
            return lines[iline], cd
        elif block == 'coord':
            istart = 8 + Latomline + Lselect
            coord = lines[istart:istart+ntotal]     # this is string of line with \n
            ### return 2D(natom by 3) coordinates
            if opt:
                #print(f'change into numbers in {whereami()}()')
                coordnum = []
                for line in coord:
                    lxyz = list(map(float, line.strip().split()))
                    coordnum.append(lxyz)
                coord=coordnum
            return coord, ntotal                    # return lines
        ### vel block: empty line followed by natoms line
        #elif block == 'vel':
        
def get_zmax(poscar):
    coords, _ = parse_poscar(poscar, block='coord', opt='number')
    #print(f'coords {coords} in {whereami()}()')
    zmax = np.array(coords).max(axis=0)[2]
    return zmax

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

def distance_pbc(p1, p2, cell_dimensions):
    """
    Calculate the Euclidean distance between two particles in a 3D periodic box (PBC).
    
    Parameters:
    p1 (np.ndarray): Coordinates of the first particle [x1, y1, z1].
    p2 (np.ndarray): Coordinates of the second particle [x2, y2, z2].
    cell_dimensions (np.ndarray): Dimensions of the periodic cell [Lx, Ly, Lz].
    
    Returns:
    float: Distance between the two particles considering periodic boundary conditions.
    """
    delta = np.array(p1) - np.array(p2)  # Vector between the two particles
    delta -= cell_dimensions * np.round(delta / cell_dimensions)  # Apply minimum image convention
    return np.linalg.norm(delta)

    # Example of periodic cell dimensions (Lx, Ly, Lz)
    #cell_dimensions = np.array([10.0, 10.0, 10.0])  # Example dimensions of the box

    # Calculate the interatomic distance considering periodic boundary conditions
    #distance_pbc = distance_pbc(particle1, particle2, cell_dimensions)

    print(f"Interatomic distance (PBC) between particle 1 and particle 2: {distance_pbc:.4f} units")
    return delta



def surf_distribution(natom, axes, cd, zmax, nlevel):
    '''For cubic axes
    natom   inserted atoms on vacuum
    zmax    max z-coord of substrate
    nlevel  distribute natoms in multiple levels
    '''
    ### how to divide the inserted atoms
    ### 2L or 3L
    dist_cret = 5

    lzoffset = []
    zoffset = 4.0               # bombing atoms to z-axis from surface, distance between O atoms
    
    Lprint = 1
    Lprintimp = 0
    ### principal axes are Ang unit
    a = np.array(axes[0])
    b = np.array(axes[1])
    c = np.array(axes[2])
    if Lprint: print(f"principal axes on x {a} y {b} z {c}")
    a_length = np.sqrt(a.dot(a))
    b_length = np.sqrt(b.dot(b))
    c_length = np.sqrt(c.dot(c))
    if Lprint: print(f"vector dot product: {a_length} {b_length} {c_length}")
    if re.match('d', cd, re.I):
        zoffset /= c_length
    for i in range(nlevel):
        lzoffset.append(zoffset*(i+1))
    print(f"reset zoffset {zoffset} due to Direct {cd}")
    zcoords = []
    for i in range(nlevel):
        print(f"{i} with {lzoffset[i]} in zoffset")
        zcoords.append(zmax + lzoffset[i])
    if Lprint: print(f"{cd}: zmax {zmax} lzoffset {lzoffset} zcoord {zcoords} in {whereami()}()")

    
    ### exclude volume for close atom
    apos = []
    bpos = []
    implant_list = []
    iatom = 0
    ilevel = 0
    i = 0
    ntotal = natom
    natom /= nlevel
    while ilevel < nlevel:
        
        while iatom < natom:
            i += 1
            if Lprintimp: print(f'{i}-th trial')
            apos = np.random.uniform(0, a_length, size=1)[0]
            bpos = np.random.uniform(0, b_length, size=1)[0]
            gen = [apos, bpos]
            if Lprintimp: print(f'{iatom+1}-th generation {gen}')
            if iatom == 0:
                implant_list.append(gen)
                iatom += 1
                if Lprintimp: print(f"implanted")
            else:
                ### compare with other atoms
                for pivot in implant_list:
                    Limplant = True
                    dist = distance_pbc(gen, pivot, [axes[0][0], axes[1][1]])     # numpy vector distance)
                    if Lprintimp: print(f"distance cret {dist_cret} < {dist} distance")
                    if dist < dist_cret:
                        Limplant = False
                        break
                if Limplant:
                    implant_list.append(gen)
                    iatom += 1
                    #print(f'generated coords {gen}')
                    if Lprintimp: print(f"implanted")
        ilevel += 1
        iatom = 0
           
    ### change coordinate to c/d
    print(f"implant list {len(implant_list)} in {whereami()}()")
    lines = []
    #print(f"cd {cd}")
    for i, xy in enumerate(implant_list):
        #lineff = ff.FortranRecordWriter('3E16.8')       #FortranRecordWriter('3E20.16')
        if re.match('d', cd, re.I):
            xy[0] /= a_length
            xy[1] /= b_length
        #print(f'x, y, z added {x} {y} {zcoord} in {whereami()}()')
        ### Now two level system
        if i < len(implant_list)/nlevel:
            line = f'{xy[0]:20.16}{xy[1]:20.16f}{zcoords[0]:20.16f}' + "\n"
        else:
            line = f'{xy[0]:20.16}{xy[1]:20.16f}{zcoords[1]:20.16f}' + "\n"
        if Lprint: print(f'formatted: {line.rstrip()}')
        lines.append(line)
    return lines

def modify_POSCAR(poscar, job='zpe', matoms=None, outf='POSCAR', option=None, nlevel=1):
    '''
    Modularize POSCAR part
    poscar      to be modified
    job     mode       s/a  select or add (job has 'add')
                zpe     s   fixedMD or selective MD for selective dynamics for ZPE calculation
                bomb    s   velocity to bombardment experiment to -z axis
                addbomb a   add molecule and velocity to -z axis
    matoms      AtomN for add atom, i (int) for selection in list
                add find lattice constand and distribute added atoms
                    append N atoms
                sel select from atoms & natoms list in POSCAR
                    ?Hf, O, Mo, S, O -> O1, O2 
                    atomlist in POSCAR for movable in zpe calculation
                velocity vel depending on T
    option      job = bomb addbomb:  T
    iatom       atom index
    atoms       atom list in POSCAR
    natoms      number of atoms list in POSCAR
    '''
    print(f"Write to {outf} in {whereami()}()")

    lines = []
    ### obtain each block and write from parse_POSCAR
    line1, atoms1 = parse_poscar(poscar, block='title'); lines.append(line1) # atoms might appear 1st line
    line, scale = parse_poscar(poscar, block='scale'); lines.append(line)
    nline, paxes = parse_poscar(poscar, block='paxes'); lines.extend(nline)
    line, atoms = parse_poscar(poscar, block='atoms')
    print(f'paxes {paxes} in {whereami()}()')
    if not line:
        line = line1
        atoms = atoms1      # line1 is saved for atom line
    ### O7 is addatom
    if matoms[0].isalpha():
        i = startnum(matoms)
        #print(f"starting number index {i} in {whereami()}() module {__name__}")
        ### 
        add_atom = matoms[:i]
        add_natom = int(matoms[i:])
    ### matoms = 3: index is select from atoms list
    else:
        ind = int(matoms) # index for selection
        
    ### add atom
    if 'add' in job:
        ### split alphabet and numeric: O7, O29, He7
        line = line.rstrip() + f'  {add_atom}' + '\n'       # do not use \t which raise type error in Vasp
        atomsold = atoms[:]
        atoms.append(add_atom)
    lines.append(line)
    line, natoms = parse_poscar(poscar, block='natoms')
    ntotalold = sum(natoms)
    if 'add' in job:
        line = line.rstrip() + f'  {add_natom}' + '\n'      # do not use \t which raise type error in Vasp
        natomsold = natoms[:]
        natoms.append(add_natom)
    lines.append(line)
    ### for sigma for MBD (Maxwell-Boltzmann distribution)
    if 'add' in job:
        atom_oldfull_list = make_atomfullist(atomsold, natomsold)       # As for original POSCAR,
    atom_fullist = make_atomfullist(atoms, natoms)
        #atom_fullist = atom_oldfull_list

    ### define number of unselected atoms (npre_unsel) and selected atoms (z-coord bombard) for selection mode
    if 'add' in job:
        npre_unsel  = ntotalold
        nselatoms   = add_natom
        ind         = -1            # for print sentence
    else:
        ### calculate index for selected atoms: movable in zpe, velocity in bombardment
        #ind = atoms.indeadd) -> not to try find using atom name
        nselatoms = natoms[ind]              # natoms to be moved for zpe
        npre_unsel = 0
        for i, na in enumerate(natoms):
            if i < ind:
                npre_unsel += na                  # npre_unsel = natoms before selected (movable) atoms
        if job == 'zpe':
            lines.append("Selective dynamics\n")
    print(f"ind {ind}  {npre_unsel} unselected in {whereami()}()")
    
    ### for Cartesian or Direct    
    line, cd = parse_poscar(poscar, block='cd'); lines.append(line)

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
        z_surf = get_zmax(poscar)
        #radius = vdw_radii[chemical_symbols.indeadd)]
        add_coords = surf_distribution(add_natom, paxes, cd, z_surf, nlevel)
        #print(f"{add_coords} in {whereami()}()")
        lines.extend(add_coords)

    ### velocity section will be provided
    ### unselected atoms have random vx, vy, vz
    ### selected atom will have -vz only with magnitue of |v|=sqrt(vx**2+vy**2+vz**2)   
    if 'bomb' in job:
        ### addtional velocity block as cartisian coordinate (A/fs)
        lines.append("\n")
        T = option
        mu = 0.0
        amukT = amu/(k_B * T)
        ms2angfs = 1.E-5         # from m/s to (1E10/1E15) Ang/fs
        lineformat = ff.FortranRecordWriter('3E16.8')
        for i, atom in enumerate(atom_fullist):
            ### mass is required: depending on mass, sigma=1/(m/k_B*T) is changed
            atomic_weight = atomic_masses[chemical_symbols.index(atom)]
            sigma = 1./np.sqrt(atomic_weight*amukT) * ms2angfs 
            vx, vy, vz = get_MBD_1D(loc=mu, scale=sigma, size=1)    # N.B. each v's are list of size 
            
            if i < npre_unsel:
                s = lineformat.write([vx[0], vy[0], vz[0]]) + "\n"
            elif i < npre_unsel + nselatoms:
                ### Define temperature for bombardment element for fast moving -> ? effective under NVT
                #amukT = amu/(k_B * 500) # T = 1000 K for fast approach without O-O dimer generation 
                #sigma = 1./np.sqrt(atomic_weight*amukT) * ms2angfs 
                #vx, vy, vz = get_MBD_1D(loc=mu, scale=sigma, size=1)    # N.B. each v's are list of size 
                v = np.sqrt(vx[0]**2 + vy[0]**2 + vz[0]**2)
                #s = lineformat.write([0.0, 0.0, -vz[0])]) + "\n"
                s = lineformat.write([0.0, 0.0, -v]) + "\n"
            else:
                s = lineformat.write([vx[0], vy[0], vz[0]]) + "\n"
            lines.append(s)
    ### print output
    filepointer = 1
    f = open (outf, 'w')
    for line in lines:
        if filepointer == 0:
            sys.stdout.write(line)
        else:
            f.write(line)
    f.close()
    return 0
