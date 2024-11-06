#!/home/joonho/anaconda3/bin/python

import argparse
import re
import ase.io.vasp as iovas

def ase_vasp_d2c(pos, job, outfile):
    '''
    Using ASE, read POSCAR/CONTCAR write POSCAR w. c/d
    pos         POSCAR
    job         d2c direct to cartesian
                c2d cartesian to direct
    outfile     POSCAR.name
    '''
    atoms = iovas.read_vasp(pos)

    if re.match('[dD]', job):
        direct=True
    else:
        direct=False

    if not outfile:
        if direct:
            outfile = pos + '_d'
        else:
            outfile = pos + '_c'
    iovas.write_vasp(outfile, atoms, direct=direct)
    
    return 0

def main():
    parser = argparse.ArgumentParser(description="add atoms, vel block")
    parser.add_argument('poscar', help="poscar to be modified")
    parser.add_argument('-j', '--job', default='d', choices=['d', 'c'], help="direct to cartesian")
    gfname =  parser.add_mutually_exclusive_group()
    gfname.add_argument('-suf', '--suffix',     help="add suffix to outfile")
    gfname.add_argument('-o', '--outfile',      help='output POSCAR name')
    args = parser.parse_args()

    ### job = bomb or addbomb
    ase_vasp_d2c(args.poscar, args.job, args.outfile)

    return 0

if __name__ == "__main__":
    main()
