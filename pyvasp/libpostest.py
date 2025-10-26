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
                opt 1, list    default(cartesian) return n * 3 2D list

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

        p_axes2D=[]   # principal axis
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
                p_axes2D.append(paxis)
                #print(f"p_axes in 2D: {paxis}")
            if opt: # return length only
                ### principal axes are Ang unit
                a = np.array(p_axes2D[0])
                b = np.array(p_axes2D[1])
                c = np.array(p_axes2D[2])
                a_length = np.sqrt(a.dot(a))
                b_length = np.sqrt(b.dot(b))
                c_length = np.sqrt(c.dot(c))
                #print(f"principle axis magnitude: {a_length} {b_length} {c_length} in {whereami()}()")
                return nline, [a_length, b_length, c_length]
            else:
                return nline, p_axes2D 
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
            ### return 2D coordinates
            if opt:
                #print(f'change into numbers in {whereami()}()')
                coordnum = []
                for line in coord:
                    lxyz = list(map(float, line.strip().split()))
                    coordnum.append(lxyz)
                coord=coordnum
            return coord, ntotal                    # return lines
        ### vel block: empty line followed by natoms line
        elif block == 'vel':
            istart = 8 + Latomline + Lselect + ntotal + 1
            try:
                b_vel = lines[istart:istart+ntotal]
            except IndexError:
                return None, None
            return b_vel, ntotal

def coord_d2c(poscar):
    '''
    in      POSCAR
    return  cartesian coordinates
    '''
    _, pvalues = parse_poscar(poscar, block='paxes', opt='length')
    coords, _ = parse_poscar(poscar, block = 'coord', opt='lis') # to 2D list
    #print(f"paxes {paxes} in function {whereami()}")
    new_coords=[]
    for coord in coords:
        acoord = []
        for i in range(len(coord)):
            xyz = coord[i] * pvalues[i]
            acoord.append(xyz)
        new_coords.append(acoord)
    return new_coords

def get_zmax(poscar):
    '''
    zmax in value
    '''
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


def obtain_atomilist0_z(zminmax, poscar, atom_species, loc):
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

def obtain_atomilist0_kind(poscar, atom_species):
    '''
    input   atom_species
    output  index list starting from 0
            return 1d, 2d
    dim     make 1d list or 2d list depending on atom species

    '''
    _, atoms  = parse_poscar(poscar, 'atoms')
    _, natoms_str = parse_poscar(poscar, 'natoms')
    natoms = list(map(int, natoms_str))

    idx1d = []
    idx2d = []
    for atom in atom_species:
        iacc = 0
        for aname, natom in zip(atoms, natoms):
            if atom == aname:
                ### start from iacc
                li = list(range(iacc, iacc+natom))
                idx1d.extend(li)
                idx2d.append(li)
            iacc += natom
    return idx1d, idx2d

    

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
    ''' obtain dirname from POSCAR[CONTCAR].name '''
    if re.match("POSCAR", poscar) or re.match("CONTCAR", poscar):
        if len(poscar) <= 7:        # just POSCAR or CONTCAR
            dirname='pos'
        ### POSCAR.dname, CONTCAR.dname
        else:
            if re.match("P", poscar):
                dirname = poscar[7:]
            elif re.match("C", poscar):
                dirname = poscar[8:]
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

def min_dist_i(p1, atoms, axes):
    dist = []
    for i, atom in enumerate(atoms):
        dist1 = distance_pbc(p1, atom, axes)
        dist.append(dist1)
        #print(f"{i}-th atom: dist = {dist1}")
        #print(f"{p1}:{atom}")
    return min(dist)

def implant_2D(pos_coords, natom, axes, cd, zfix, zmax, r_crit=3.0, nlevel=1):
    '''For cubic axes
    pos_coords  original coords of POSCAR in cartesian for pbc comparison
    natom   inserted atoms on vacuum
    zfix  list with elements: 'top' make a distance
                                one z-values for fixed position
                                two z-values for inbetween
    zmax    max in value (C|D)
    r_crit  implantation criteria for atomic inter-distance
            donot apply for added atoms in case smaller value for planar implant
    nlevel  distribute natoms in multiple levels
    '''
    ### how to divide the inserted atoms
    ### 2L or 3L
    if re.match('t', zfix[0]):
        ztag = 'top'
    else:
        ztag = 'inter'
        zmax = 0
    lzoffset = []
    if ztag == 'top':
        zoffset = 4.0               # bombing atoms to z-axis from surface, distance between O atoms
        interdist = 5               # compare with inserted O
    else:
        zoffset = 0.0
        interdist = r_crit               # compare with all atoms
        interdist_addatom = 4
    
    print(f"ztag = {ztag}")
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

    ### Convert to C
    if re.match('d', cd, re.I):
        #zoffset /= c_length             # converted to D
        zmax *= c_length                # converted to C
    for i in range(nlevel):
        lzoffset.append(zoffset+(zoffset+1)*i)
    print(f"reset zoffset {zoffset} due to Direct {cd} in function {whereami()}()")


    ### zcoordinates: use zmax for top
    ### convert to cartesian from direct
    if ztag == 'top':
        zcoord = zmax                   # this is C
    ### use zfix for z-coordinates
    else:
        if len(zfix) == 1:
            zcoord = float(zfix[0])
        else:
            zcoord = (float(zfix[0]) + float(zfix[1]))/2.     # inbetween
    print(f"zcoord {zcoord}")
    zcoords = []
    ### for inter model no nlevel
    for i in range(nlevel):
        print(f"{i} with {lzoffset[i]} in zoffset")
        if ztag == 'top':
            zcoords.append(zcoord + lzoffset[i])
        else:
            zcoords.append(zcoord)
    if Lprint: 
        print(f"{cd}: zmax {zmax} lzoffset {lzoffset} zcoord {zcoords} in {whereami()}()")
        print(f"if takes longer time, reduce interatomic distance by -d less than {r_crit}")

    if ztag == 'top':
        comp_atoms = []             # no need to compare with existing atoms
    else:
        comp_atoms = pos_coords     # compare with existing atoms

    ### exclude volume for close atom
    apos = []
    bpos = []
    implant_list = []
    i_try = 0
    ntotal = natom
    natom_level = natom/nlevel                     # natom in a level
    iatom = 0
    ilevel = 0

    ### ilevel loop in case added atoms in multilevel in z-coord
    while ilevel < nlevel:
        #if ztag == 'inter':
        #    implant_list.extend(coord)
        ### calculation using Cartesian
        natom_orig = len(pos_coords)
        added_atom_coords=[]
        ### natom loop
        while iatom < natom_level:
            i_try += 1
            if Lprintimp: print(f'{i}-th trial')
            apos = np.random.uniform(0, a_length, size=1)[0]
            bpos = np.random.uniform(0, b_length, size=1)[0]
            gen = [apos, bpos, zcoords[ilevel]]          # generation in cartesian
            if Lprintimp: print(f'{iatom+1}-th generation {gen}')
            ### check distance
            if iatom == 0 and ztag == 'top' and ilevel == 0:
                comp_atoms.append(gen)
                iatom += 1
                if Lprintimp: print(f"implanted")
            ### compare with O's for top in 2D, all for fixed in 3D
            else:
                ### compare with other atoms
                ### use different crit for other atoms and adding atoms
                for i, pivot in enumerate(comp_atoms):
                    Limplant = True
                    dist = distance_pbc (gen, pivot, [axes[0][0], axes[1][1], axes[2][2]])
                    #if Lprintimp: print(f"distance cret {interdist} < {dist} distance")
                    if ztag == 'top':
                        if dist < interdist:
                            Limplant = False
                            break
                    else:
                        if i < natom_orig:
                            dist_cret = interdist
                        else:
                            dist_cret = interdist_addatom
                        if dist < dist_cret:
                            Limplant = False
                            break
                if Limplant:
                    ### check min distance
                    #print(f"min dist: {min_dist_i(gen, comp_atoms,[axes[0][0], axes[1][1], axes[2][2]])}")
                    comp_atoms.append(gen)
                    #added_atom_coords.append(gen)
                    iatom += 1
                    #print(f'generated coords {gen} in loop: iatom < natom_level')
                    if Lprintimp: print(f"implanted")
        iatom = 0
        ilevel += 1
           
    ### change coordinate to c/d
    print(f"implant list {len(implant_list)-natom_orig} in {whereami()}()")
    lines = []
    print(f"cd {cd}, total natom {ntotal}")
    #print(f"{comp_atoms[-ntotal:]}")
    #print(f"ntotal: {ntotal} size {len(comp_atoms)}")
    for i, xy in enumerate(comp_atoms[-ntotal:]):
        lineff = ff.FortranRecordWriter('3E16.8')       #FortranRecordWriter('3E20.16')
        if re.match('d', cd, re.I):
            xy[0] /= a_length
            xy[1] /= b_length
            xy[2] /= c_length
        #print(f'x, y, z added {x} {y} {zcoord} in {whereami()}()')
        ### Now two level system
        #if i < len(implant_list)/nlevel:
        #if i < len(comp_atoms) - natom_orig:
        line = f'{xy[0]:20.16}{xy[1]:20.16f}{xy[2]:20.16f}' + "\n"
        #else:
        #    line = f'{xy[0]:20.16}{xy[1]:20.16f}{zcoords[1]:20.16f}' + "\n"
        if Lprint: print(f'formatted: {line.rstrip()}')
        lines.append(line)
    return lines

def modify_POSCAR(poscar, job='zpe', mode=None, mod_atoms=None, zpos=None, temp=300, htemp=None,\
            vel_type='random',v_reverse=False, outf='POSCAR', r_crit=None, asort=None, nlevel=None):
    '''
    Modularize POSCAR part
    inputs:
        poscar      to be modified
        job         zpe     s(i)   fixedMD or selective MD for selective dynamics for ZPE calculation
                    bomb    s(l)   velocity to bombardment experiment to -z axis
                    addbomb a   add molecule and velocity to -z axis
                    vel      N/A assign velocity in the given configuration
                    sort    sl    change latomnames, lnatoms
                    rm      si   atom index
        mode        sl, si  selection index in atom list or just atom index
                    a       add format O12; atom symbol:natom
        mod_atoms   modified atoms
                    sl       atoms in the list index
                    a       add find lattice constand and distribute added atoms
                            append N atoms
                    s   select from atoms & natoms list in POSCAR
                        ?Hf, O, Mo, S, O -> O1, O2 
                        atomlist in POSCAR for movable in zpe calculation
                    vel assign velocity to all the atoms depending on T
        zpos        the posotion in z-axis where atoms to be added
                    top: above of surface atom 4 A away from top atom
                    z1 [z2]: position of z-value or inbetween two z-values
        temp        temperature for velocity  T
        htemp       hyper_temperature for hyperthermal species in plasma
        vel_type    'random' for assign following T
                    'copy' to copy original file
        v_reverse   boolean for bombing to +z direction
        nlevel      number of levels to add O atoms in multi level
        r_crit      reference distance between atoms: as for existing atoms not generating atoms
    variables:
        iatom       atom index
        atoms       atom list in POSCAR
        natoms      number of atoms list in POSCAR
        matoms      a   atom species followed by natom to be added in modifying POSCAR
                    sl   integer in atom species line
                        two integers linked by '-' to be contract in atom list
                    si   integers of atom list
    '''
    ### define constants
    #nlevel = 1         # deprecated for now
    z_top = 4           # Ang from top surface
    interdist = 4       # between added O atoms
    interdist_mid = 3   # at crowded space of interface

    print(f"Write to {outf} in {whereami()}()")
    ### obtain (modify) atom indices
    print(f"mod_atoms {mod_atoms} {len(mod_atoms)}")
    if len(mod_atoms) == 1:
        mod_atoms0 = mod_atoms[0]   # use mod_atoms0 for 'sl' or 'a'

    lines = []
    ### obtain each block and write from parse_POSCAR
    line1, atoms1 = parse_poscar(poscar, block='title'); lines.append(line1) # atoms might appear 1st line
    line, scale = parse_poscar(poscar, block='scale'); lines.append(line)
    nline, paxes = parse_poscar(poscar, block='paxes'); lines.extend(nline)

    ### deal with atomname & natom lines together to deal with job=rm
    line, latoms = parse_poscar(poscar, block='atoms')
    print(f'paxes {paxes} in {whereami()}()')
    line_na, lnatoms = parse_poscar(poscar, block='natoms')
    ntotalold = sum(lnatoms)
    
    print(f"mode: {mode}, job: {job} in {whereami()}()")
    if not line:
        line = line1
        latoms = atoms1      # line1 is saved for atom line
    ### O7 is addatom
    if mode == 'a': #mod_atoms0.isalpha():
        i = startnum(mod_atoms0)
        #print(f"starting number index {i} in {whereami()}() module {__name__}")
        add_atomname = mod_atoms0[:i]           # atom name to be added
        add_natom = int(mod_atoms0[i:])     # natom to be added
        ### split alphabet and numeric: O7, O29, He7
        line = line.rstrip() + f'  {add_atomname}' + '\n'       # do not use \t which raise type error in Vasp
        latomsold = latoms[:]
        latoms.append(add_atomname)
        ### deal with natom line
        line_na = line_na.rstrip() + f'  {add_natom}' + '\n'      # do not use \t which raise type error in Vasp
        natomsold = lnatoms[:]
        lnatoms.append(add_natom)

    ### matoms is integer in natom lines of POSCAR for mode='sl' or atom index list for mode='si'
    elif mode == 'sl':
        if mod_atoms0.isdigit():
            ind = int(mod_atoms0)               # index for atom selection in POSCAR atom line
        elif job == 'sort' and re.search('-',matoms):
            indices = re.split('-', matoms)
            ind = int(indices[0])
            indf = int (indices[-1])
            ### contract atom line
            atom_kinds_tobesorted = latoms[ind:indf+1]
            if asort:
                atom_sort = asort
            else:
                print(f"input atoms to be sorted: -as ")
                atom_sort = get_atom_kinds(atom_kinds_tobesorted)
            #print(f"atom_sort {atom_sort}")
            latom_indices_tobesorted=[ atom_sort.index(a) for a in atom_kinds_tobesorted]
            new_latoms = latoms[:ind]
            new_latoms.extend(atom_sort)
            new_latoms.extend(latoms[indf+1:])
            print(f"{new_latoms}")
            line = "  ".join(new_latoms) + "\n"
            #print(line); print(f"latom line {latom_indices_tobesorted}")
            #if job == 'sort':               # change natoms list
            ### sum natomslist following latom_indices_tobesorted
            lnatom_tobesorted = lnatoms[ind: indf+1]
            nselatoms = sum(lnatom_tobesorted)
            new_natoms_tobesored = np.zeros(len(atom_sort)).astype(int)
            for i, inatom in enumerate(lnatom_tobesorted):
                #print(f"{i}: {inatom}")
                new_natoms_tobesored[latom_indices_tobesorted[i]] += inatom
            new_natoms = lnatoms[:ind]
            new_natoms.extend(list(new_natoms_tobesored))
            new_natoms.extend(lnatoms[indf+1:])
            #print(f"new natoms: {new_natoms}")
            line_na = "  ".join(map(str, new_natoms)) + "\n"
            #print(f"natom line: {new_natoms_tobesored}"); # sys.exit(12)
            #print(f"lnatom_tobesorted: {lnatom_tobesorted}")
    elif mode == 'si': #selection by index
        ### change str into integer
        int_matoms = list(map(int, mod_atoms))
        if job == 'rm':
            rm_lnatoms_ind=[]
            for iatom in int_matoms: # iatom is sorted
                sum_natoms = 0
                for i, natom in enumerate(lnatoms):
                    sum_natoms += int(natom)
                    if iatom < sum_natoms:   # delete atom
                        rm_lnatoms_ind.append(i)
                        break
            for ind_rm in rm_lnatoms_ind:
                lnatoms[ind_rm] -= 1
                if lnatoms[ind_rm] == 0:
                    print(f"all the {ind_rm}-th atoms are removed")
                    sys.exit(100)
            line_na = "  ".join(map(str, lnatoms)) + "\n"
            print(f"line_na {line_na} in {whereami()}()")
    elif mode == 'vel':
        ind = -1                        # all the atoms in atom index
    ### add atom
    lines.append(line)                  # job=rm does not remove all the atoms
    #sys.exit(11)
    ### deal with natom line
    lines.append(line_na)
    print(line_na); 
    
    # ind is for natoms line index
    #print(f"ind {ind}  {npre_unsel} unselected in {whereami()}()")

    ### for sigma for MBD (Maxwell-Boltzmann distribution)
    #if 'add' in job:
    #    atom_oldfull_list = make_atomfullist(atomsold, natomsold)       # As for original POSCAR,
    atom_fullist = make_atomfullist(latoms, lnatoms)
        #atom_fullist = atom_oldfull_list

    ### define number of unselected atoms (npre_unsel) and selected atoms (z-coord bombard) for selection mode
    if mode == 'a':
        npre_unsel  = ntotalold
        nselatoms   = add_natom
        ind         = -1            # for print sentence
    elif mode == 'sl':
        ### calculate index for selected atoms: movable in zpe, velocity in bombardment
        #ind = atoms.indeadd) -> not to try find using atom name
        nselatoms = lnatoms[ind]              # natoms to be moved for zpe
        npre_unsel = 0
        for i, na in enumerate(lnatoms):
            if i < ind:
                npre_unsel += na                  # npre_unsel = natoms before selected (movable) atoms
        if job == 'zpe':
            lines.append("Selective dynamics\n")
    elif mode == 'vel':
        npre_unsel  = ntotalold
        nselatoms = 0
        

    
    ### for Cartesian or Direct    
    line, cd = parse_poscar(poscar, block='cd'); lines.append(line)

    ### Deal with COORDINATEs in line list
    coords, _ = parse_poscar(poscar, block='coord') 
    
    ## modify selective dynamics for ZPE
    new_coords=[]
    new_coords_rem=[]
    if job == 'zpe':
        ### coords -> new_coords in ZPE
        for i, line_coord in enumerate(coords):
            if i < npre_unsel:
                new_line = line_coord.rstrip() + " F F F\n"
            elif i < npre_unsel + nselatoms:
                new_line = line_coord.rstrip() + " T T T\n"
            elif i < ntotalold:
                new_line = line_coord.rstrip() + " F F F\n"
            new_coords.append(new_line)
        coords = new_coords
    elif job == 'sort':
        ### directly divide coords by indices
        new_sort_2Dcoords = [ [ ] for i in range(len(atom_sort))]
        #for i in range(len(atom_sort)):
        #    new_coords_sort[i] = []
        jsort = 0
        iatom = 0
        for i, natom in enumerate(lnatoms):
            if i < ind:
                new_coords.extend(coords[iatom:iatom+natom])
            elif i < indf:
                new_sort_2Dcoords[latom_indices_tobesorted[i-ind]].extend(coords[iatom:iatom+natom])
            else:
                ### treat sorted coordinates and add the remainder
                new_coords_rem.extend(coords[iatom:iatom+natom])
            iatom += natom
        ### combind coords here
        for i in range(len(atom_sort)):
            new_coords.extend(new_sort_2Dcoords[i])
        new_coords.extend(new_coords_rem)
        lines.extend(new_coords)
    elif job == 'rm':
        print(f"len coords {len(coords)} in {whereami()}()")
        for i in int_matoms:
            coords = coords[:i] + coords[i+1:]
        print(f"len coords {len(coords)} in {whereami()}()")
        lines.extend(coords)
    ### job == 'bomb' or 'add'
    else:
        lines.extend(coords)                            # copy original coords
        ### Make additional coordinates for added atoms
        if mode == 'a':
            ### for Orthorhombic crystal structure
            #z_coord = 'above'
            zmax = get_zmax(poscar)
            #radius = vdw_radii[chemical_symbols.indeadd)]
            #if zpos:
            #    z_line = zpos
            #    z_coord = 'fixed'
            if re.match('d', cd, re.I):
                d2coords_cart = coord_d2c(poscar)
            else:
                d2coords_cart, _ = parse_poscar(poscar, block='coord', opt='lis')

            print(f"{d2coords_cart[0]} in function {whereami()}()") 
            add_coords = implant_2D(d2coords_cart, add_natom, paxes, cd, zpos, zmax, r_crit, nlevel)
            #print(f"{add_coords} in {whereami()}()")
            lines.extend(add_coords)

    ### velocity section will be provided
    ### unselected atoms have random vx, vy, vz
    ### selected atom will have -vz only with magnitue of |v|=sqrt(vx**2+vy**2+vz**2)   
    ### for job == 'rm' use 'rmmd' to include vel
    if 'bomb' in job or 'md' in job or vel_type:
        vel_orig, _ = parse_poscar(poscar, block='vel')
        ### addtional velocity block as cartisian coordinate (A/fs)
        print(f"original vel {len(vel_orig)} in job {job}")
        if job == 'rm' and vel_orig:
            for i in int_matoms:
                vel_orig = vel_orig[:i] + vel_orig[i+1:]
            print(f"len velocity {len(vel_orig)} in {whereami()}()")
            lines.append(f"{os.linesep}")
            lines.extend(vel_orig)
        else:
            print(f"Add velocity at {temp} K for up to {npre_unsel} atom")
            lines.append("\n")
            T = temp
            mu = 0.0                # mean value for M-B distribution
            amukT = amu/(k_B * T)
            ms2angfs = 1.E-5         # from m/s to (1E10/1E15) Ang/fs
            lineformat = ff.FortranRecordWriter('3E16.8')
            for i, atom in enumerate(atom_fullist):
                ### mass is required: depending on mass, sigma=1/(m/k_B*T) is changed
                atomic_weight = atomic_masses[chemical_symbols.index(atom)]
                sigma = 1./np.sqrt(atomic_weight*amukT) * ms2angfs 
                vx, vy, vz = get_MBD_1D(loc=mu, scale=sigma, size=1)    # N.B. each v's are list of size 
                ### if POSCAR has velocity block, read it
                if i < npre_unsel:
                    if re.match('c', vel_type):
                        s = vel_orig[i]
                    else:
                        s = lineformat.write([vx[0], vy[0], vz[0]]) + "\n"
                elif i < npre_unsel + nselatoms:
                    ### Define temperature for bombardment element for fast moving -> ? effective under NVT
                    #amukT = amu/(k_B * 500) # T = 1000 K for fast approach without O-O dimer generation 
                    #sigma = 1./np.sqrt(atomic_weight*amukT) * ms2angfs 
                    #vx, vy, vz = get_MBD_1D(loc=mu, scale=sigma, size=1)    # N.B. each v's are list of size 
                    ### in case hyperthermal species, convert eV to T and assign to selected atoms
                    print(f"{vx[0]:6.3f} {vy[0]:6.3f} {vz[0]:6.3f}")

                    if htemp and htemp < 100.:
                        # use eV in velocity units
                        pass
                    elif htemp: # use high T as velocity units
                        T = htemp
                        amukT = amu/(k_B * T)
                        sigma = 1./np.sqrt(atomic_weight*amukT) * ms2angfs
                        vx, vy, vz = get_MBD_1D(loc=mu, scale=sigma, size=1)
                    #else: # use the same T as substrate
                    v = np.sqrt(vx[0]**2 + vy[0]**2 + vz[0]**2)
                    #s = lineformat.write([0.0, 0.0, -vz[0])]) + "\n"
                    if 'bomb' in job:
                        if not v_reverse:
                            v *= -1
                        s = lineformat.write([0.0, 0.0, v]) + "\n"
                    else:
                        s = lineformat.write([vx[0], vy[0], vz[0]]) + "\n"
                else:
                    if re.match('c', vel_type):
                        s = vel_orig[i]
                    else:
                        s = lineformat.write([vx[0], vy[0], vz[0]]) + "\n"
                #print(f"{s} in function {whereami()}()")
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
