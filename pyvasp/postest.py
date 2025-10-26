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
from libpostest import modify_POSCAR

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
    job bomb    make velocity with v_z for bombing atoms
        add     increase atoms
        zpe
        md      just include velocity without bombing
        sort    sort atoms, natoms line by two indices
    #Depending on Job modify these arguments
    lnatom  add_atoms|sel_list
        a[dd]   increas the natom
                'O4' atom name followed by number of atoms to be added
            (bomb, add)     
        s[el]   select atoms by index in atoms line
        e.g.:
            (bomb)  e.g. 3      : one index which starts from 0
            (sort)  e.g. 1-10   : two indices linked by dash
    lmd     vel     include velocity section 
    vtype  [random(default)| -z | +z]
            copy from input POSCAR
    t[emp]      system temperature
    ht[emp]     hyper temperature for inserted atoms
            top: increase simulation time
            interface: overcome attraction from both side
                to determine the interacting surface by initial v_z
    '''
    parser = argparse.ArgumentParser(description="add atoms, vel block")
    parser.add_argument('poscar', default='POSCAR', help="poscar to be modified")
    parser.add_argument('-j', '--job', default='add', choices=['bomb','add','rm','zpe', 'md','sort'],\
                        help="VASP job for POSCAR")
    ### depeselect existing atom or add atoms for bomb
    #gatoms =  parser.add_mutually_exclusive_group()
    parser.add_argument('-m', '--mode', choices=['sl', 'si', 'a', 'vel'], help="mode: selection from atom list index or atom index or add atoms")
    parser.add_argument('-a', '--atoms', nargs='*', help="list index, atom index or natom")
    parser.add_argument('-z', '--zcoord', default = ['top'], nargs='*', help="'top', one or two z-coord")
    parser.add_argument('-d', '--distance', default = 3.0, type=float, help="interdistance creteria for implantation")
    parser.add_argument('-t', '--temp', type=float, default=500,  help="T(K) for atomic velocity")
    parser.add_argument('-ht', '--hypertherm', type=float, help="T(eV), for atom velocity")
    parser.add_argument('-vt', '--velocity_type', default='random', choices=['random', 'copy'], help="T for atom velocity, why not random")
    parser.add_argument('-vr', '--vel_reverse', action='store_true', help="make bombing to upside")
    parser.add_argument('-l', '--nlevel', type=int, default=1,  help="atoms displaced in multi levels")
    parser.add_argument('-as', '--sort', nargs='*', help="order of atoms in sorting")
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
                \n\t\tkpy pos_modify.py CONTCAR.HfSe2L1 -j bomb -a O8 -t 800 -ht 1000 -o POSCAR.HfSe2L1 -z 10 -d 2.5\
                \n\t    rm  : mode = 'si'\
                \n\t-m mode:\
                \n\t    sl: -a 3-9, select list interval\
                \n\t    si: -a 3-7 or 1 3 8 select atom indices\
                \n\t    a : -a O6 atom_name + natom\
                ")
        sys.exit(1)
        # if no selection, assign vel to all atoms in given config
    if args.mode:
        mode = args.mode
    else:
        if args.job == 'sort':
            mode = 'sl'     # select atom list interval: 3-7
        elif args.job == 'rm':
            mode = 'si'     # select atom index: 3-7 or 3 5 8 etc
        elif args.job == 'bomb' or args.job == 'add': # line for example
            mode = 'a'      # add append 
        else:
            mode = 'vel'
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
        outfile = args.poscar + args.job

    #if not args.velocity_type:
    #    if 'md' in args.job:
    #        vel_type = 'random'
    #    elif 'bomb' in args.job or args.job == 'rm':
    #        vel_type = 'copy'
    #if 'bomb' in args.job:
    ### job = bomb or addbomb
    #pos_bombardment(args.poscar, args.job, atoms, args.zcoord, args.temp, args.velocity, args.nlevel, outfile)
    modify_POSCAR(args.poscar, job=args.job, mode = mode, mod_atoms=args.atoms, zpos=args.zcoord, \
    temp=args.temp, htemp=args.hypertherm, vel_type=args.velocity_type, v_reverse=args.vel_reverse, nlevel=args.nlevel,\
    asort=args.sort, r_crit=args.distance, outf=outfile)

    return 0

if __name__ == "__main__":
    main()
