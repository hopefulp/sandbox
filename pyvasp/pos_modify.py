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
    job bomb    make velocity with v_z for bombing atoms
        add     increase atoms
        zpe
        md      just include velocity without bombing
        sort    sort atoms, natoms line by two indices
    add_atoms|sel_atoms
        a       'O4' atom name followed by number of atoms to be added
            (bomb, add)     
        s       select atoms by index in atoms line
            (bomb)  e.g. 3      : one index which starts from 0
            (sort)  e.g. 1-10   : two indices linked by dash
    velocity_t  (default) random by -t and -ht
                copy from input POSCAR
    '''
    parser = argparse.ArgumentParser(description="add atoms, vel block")
    parser.add_argument('poscar', help="poscar to be modified")
    parser.add_argument('-j', '--job', default='add', choices=['bomb','add','zpe', 'md','sort'], help="job of poscar changing")
    ### select existing atom or add atoms for bomb
    gatoms =  parser.add_mutually_exclusive_group()
    gatoms.add_argument('-s', '--sel_atoms', help="one atom species or index in POSCAR: Hf O Mo S O  0 1 2 not yet for 2-5")
    gatoms.add_argument('-a', '--add_atoms', help="add atoms: O12 Fe3 3-index")
    parser.add_argument('-z', '--zcoord', default = ['top'], nargs='*', help="'top', one or two z-coord")
    parser.add_argument('-d', '--distance', default = 3.0, type=float, help="interdistance creteria for implantation")
    parser.add_argument('-t', '--temp', type=float, default=500,  help="T(K) for atomic velocity")
    parser.add_argument('-ht', '--hypertherm', type=float, help="T(eV), for atom velocity")
    parser.add_argument('-vt', '--velocity_type', default='random', choices=['random', 'copy'], help="T for atom velocity")
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
                \n\t    kpy pos_modify.py POSCAR.HfSe2L1 -j sort -s 1-17 -as Se O\
                \n\tbomb:: input velocity for -a addatoms or -s select by index with -t temp and -ht hypertemp\
                \n\t    kpy pos_modify.py CONTCAR.HfSe2L1 -j bomb -a O8 -t 800 -ht 1000 -o POSCAR.HfSe2L1 -z 10 -d 2.5\
                ")
        sys.exit(1)
    if args.sel_atoms:
        atoms='s' + args.sel_atoms
    elif args.add_atoms:
        atoms='a' + args.add_atoms
    else:
        # if no selection, assign vel to all atoms in given config
        atoms='vel'
        #print(f"input error: select or add atoms via -s or -a")
        #sys.exit(1)
    ### outfile name need to be passed
    if args.outfile:
        outfile = args.outfile
    elif args.suffix:
        outfile = args.poscar + args.suffix
    else:
        outfile = args.poscar + args.job

    #if 'bomb' in args.job:
    ### job = bomb or addbomb
    #pos_bombardment(args.poscar, args.job, atoms, args.zcoord, args.temp, args.velocity, args.nlevel, outfile)
    modify_POSCAR(args.poscar, job=args.job, mode_atoms=atoms, zpos=args.zcoord, \
    temp=args.temp, htemp=args.hypertherm, vel_type=args.velocity_type, v_reverse=args.vel_reverse, nlevel=args.nlevel,\
    asort=args.sort, r_crit=args.distance, outf=outfile)

    return 0

if __name__ == "__main__":
    main()
