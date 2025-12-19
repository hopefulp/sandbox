#!/home/joonho/anaconda3/bin/python
'''
add atom
add velocity block
'''
import argparse
import sys
#import re

#from common import whereami
#import chem_space as cs
from libposcar import modify_POSCAR

def pos_bombardment(pos, job, atoms, zposition, temp, vel, nlevel, outfile):
    '''
    pos         POSCAR
    job         add   -> just append atoms to POSCAR
                bomb  -> atoms have velocity
                
    atoms   a.. add atoms
                aO20    add 20 O's
            s..  select atoms in POSCAR
                sO2     select 2nd O group in atom list of POSCAR
                sOn     select 1st O group and n atoms
    outfile     POSCAR.name
    '''
    
    #if re.match('s', bomb_atoms):
    modify_POSCAR(pos, job=job, mode_atoms=atoms, zpos=zposition, temp=temp, outf=outfile)
    #elif re.match('a', bomb_atoms):
    #    add_atoms(pos, bomb_atoms[1:])
    
    return 0

def main():
    '''
    job bomb    make velocity with v_z for bombing atoms: -
        add     increase atoms
        zpe
        md      just include velocity without bombing
        sort    sort atoms, natoms line by two indices
    #Depending on Job modify these arguments
    natom       add_atoms|sel_list
    Latom: 'a' add or 's' select
        a[dd]   increas the natom
                'O4' atom name followed by number of atoms to be added
            (bomb, add)     
        s[el]   select atoms by index in atoms line
        e.g.:
            (bomb)  e.g. 3      : one index which starts from 0
            (sort)  e.g. 1-10   : two indices linked by dash
    lmd     vel     include velocity section 
    vtype  [r[andom]|copy|-z|+z]
            copy from input POSCAR
    t[emp]      system temperature
    ht[emp]     hyper temperature for inserted atoms
            top: increase simulation time
            interface: overcome attraction from both side
                to determine the interacting surface by initial v_z
    '''
    parser = argparse.ArgumentParser(description="add atoms, vel block")
    parser.add_argument('poscar', default='POSCAR', help="poscar to be modified")
    parser.add_argument('-j', '--job', default='add', choices=['bomb','inter','add','rm','zpe','md','sort','split','atype'],\
                        help="VASP job for POSCAR")
    ### select and add might be compatible
    gatoms =  parser.add_mutually_exclusive_group()
    gatoms.add_argument('-s', '--select', help="[l,i,aname]; l:atomname list, i:atomindex list, atomname: l6,l6-4,l6+3,O,O-4,i3-7 etc")
    #parser.add_argument('-sn', '--aselect_number', help="if select atom atom kind, number of atoms")
    gatoms.add_argument('-a', '--add', help="atom kinds followed by natom O4 S3 etc")
    #gatoms.add_argument('-a', '--add', nargs='*', help="atom kinds followed by natom O4 S3 etc")
    gvel = parser.add_argument_group(title='velocity')
    gvel.add_argument('-v', '--Lvelocity', action='store_true', help="include velocity in POSCAR")
    gvel.add_argument('-t','-T', '--temp', type=float, default=500,  help="system temperature T(K) for atomic velocity")
    gvel.add_argument('-ht', '--hypertherm', default=600, type=float, help="velocity for input atoms T(eV,K)")
    gvel.add_argument('-vt', '--vel_type', default='zdn', choices=['r','random', 'copy', 'zup', 'zdn'], help="T for atom velocity, why not random")
    gvel.add_argument('-vr', '--vel_reverse', action='store_true', help="make bombing to upside")
    parser.add_argument('-z', '--zcoord', default = ['top'], nargs='*', help="'top', one or two z-coord")
    parser.add_argument('-d', '--distance', default = 3.0, type=float, help="interdistance creteria for implantation")
#    parser.add_argument('-l', '--nlevel', type=int, default=1,  help="atoms displaced in multi levels")
    ### may combine job-dependent unique option:: sort - sort list, 
    parser.add_argument('-sl', '--sort_list', nargs='*', help="sort list: order of atoms in sorting")
    parser.add_argument('-u', '--usage', action='store_true', help = 'print usage')
    gfname =  parser.add_mutually_exclusive_group()
    gfname.add_argument('-suf', '--suffix',     help="add suffix to outfile")
    gfname.add_argument('-o', '--outfile',      help='output POSCAR name')
    args = parser.parse_args()

    if args.usage:
        print(f"\tsort:: by atom line index interval\
                \n\t    kpy pos_modify.py POSCAR.HfSe2L1 -j sort -a 1-17 -as Se O\
                \n\t-j job = bomb, add, rm, zpe, md, sort\
                \n\t    sort: mode = 'sl', -a list interval\
                \n\t    bomb: mode = 'sl' or 'a', input velocity for -a addatoms or -a select by index with -t temp and -ht hypertemp\
                \n\t    split: aselect index lO\
                \n\t\tkpy pos_modify.py CONTCAR.HfSe2L1 -j bomb -a O8 -t 800 -ht 1000 -o POSCAR.HfSe2L1 -z 10 -d 2.5\
                \n\t    rm  : mode = 'si'\
                \n\t    -s:  atom selection: O, O-4, O+4, l2, i3-9, select list interval\
                \n\t\tO: oxygen, O-4: split O from last 4, l1: second atom, i3-7: atom index\
                \n\t    -a : add atom, -a O6 atom_name + natom\
                ")
        sys.exit(1)
        # if no selection, assign vel to all atoms in given config
    ### args.mode is divided into -s and -a
    ### -s and -a is not compatible at the moment

    ### if atoms has 'S': select in POSCAR else add to POSCAR
    if args.add:
        atoms = 'A' + args.add
    else:
        if args.select:
            atoms = args.select
        else:
            atoms = 'all'       # if not add, not select, select all the atoms

    ### Using job: assign 'S' or 'A'
    Lvelocity = args.Lvelocity
    if args.job == 'sort' or args.job == 'md' or args.job == 'rm':
        if args.job == 'md':
            Lvelocity = True
    elif args.job == 'bomb' or args.job == 'add' or args.job == 'inter': # line for example
        if args.job == 'bomb' or args.job == 'inter':
            Lvelocity = True
            if args.job == 'bomb':
                vel_type = 'zdn'
            elif args.job == 'inter':
                vel_type = 'random'
    if args.vel_type:
        vel_type = args.vel_type

    ### make temp dict
    dic_vel = {}
    if Lvelocity:
        dic_vel={'t':args.temp, 'ht': args.hypertherm, 'vt': args.vel_type,'vr':args.vel_reverse}
    #print(f"input error: select or add atoms via -s or -a")
    #sys.exit(1)
    ### outfile name need to be passed
    if args.outfile:
        outfile = args.outfile
        if not 'POSCAR' in outfile:
            outfile = 'POSCAR.' + outfile
    elif args.suffix:
        outfile = args.poscar + args.suffix
    else:
        if args.job == 'split' or args.job == 'atype':
            outfile = args.poscar + ".vas"
        else:
            outfile = args.poscar + args.job

    ### nlevel=args.nlevel,
    modify_POSCAR(args.poscar, job=args.job, xatoms=atoms, dic_vel=dic_vel, zpos=args.zcoord,\
    asort=args.sort_list, r_crit=args.distance, outf=outfile)

    return 0

if __name__ == "__main__":
    main()
