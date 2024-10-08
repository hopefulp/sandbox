#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os
import numpy as np
from _extract_line import *
from common import *
from ase import Atoms, Atom
import ase.io

### scratch keywords for AMP from IM/MM input file of the form of .fin

kw_natom = "natm"
kw_npts = "npts"
kw_ene = "energy"
kw_coord = "cartesian"

def extract_data(f_input, f_atom, fit_e):

    fname = fname_parsing(f_input)
    fname_traj = fname[0] + '.traj'
    traj = ase.io.Trajectory(fname_traj, 'w')
    atom_list=[]
    ##atoms = Atoms()
    fout=open("ene.dat", 'w')
    tag_ene = "OFF"
    tag_coord = "OFF"

    natom = int(extract_one_line(f_input, kw_natom).rstrip())
    npoint = int(extract_one_line(f_input, kw_npts).rstrip())
    print("natom = %d, npoints = %d" % (natom, npoint))
    atom_name = extract_col_qchem(f_atom, "molecule")
    print(atom_name, len(atom_name))
    energies = []
    coordinates = []

    i=0
    with open(f_input, 'r') as f:
        for line in f:          # line has "\n"
            #print line,         # print writes its own "\n"
            if tag_ene == "OFF":
                if re.search(kw_ene, line):
                    tag_ene = "ON"
                    tag_write_ene = "OFF"
                    i+=1
            # inside energy tag                    
            elif tag_ene == "ON":
                #print line,
                # obtain energy -> 1 line
                if tag_write_ene == "OFF":
                    tag_write_ene = "ON"
                    atoms = Atoms()                             # initialize Atoms object
                    energy = float(line.rstrip())
                    energies.append(energy)
                    fout.write(line)    # write does not write "\n"
                    #print "%d-th energy: %s" % (i, line.rstrip())
                if tag_coord == "OFF":
                    if re.search(kw_coord, line):
                        mol_coord=[]
                        tag_coord = "ON"
                        j = 0
                elif tag_coord == "ON":
                    j += 1              # read next coordinates line
                    lcoord = line.split()
                    mol_coord.append(lcoord)
                    ai = Atom(atom_name[j-1], lcoord)           # Make Atom object
                    atoms.append(ai)                            # make Atoms object
                    #print lcoord
                    if j == natom:
                        tag_coord = "OFF"
                        tag_ene = "OFF"
                        coordinates.append(mol_coord)
                        atoms.set_potential_energies(energy)
    fout.close()
    #print coordinates
    #print "%d %d" % (len(energies), len(coordinates))
    if i != npoint or i != len(energies) or i != len(coordinates):
        print("scratch %d geometries are failed" % i)
    else:
        print("%d geometries are scratched" % i)

    ### np.save
    np_e = np.array(energies)
    np_c = np.array(coordinates)
    np.savez("np_e_c", ene=np_e, coord=np_c)
    return

def main():
    parser = argparse.ArgumentParser(description='get coordinate [energy, force, hessian] from IM input file (.fin)')
    parser.add_argument('f_immm_input', help='IM/MM input file')
    parser.add_argument('-a', '--atom_file', help='Atom list')
    parser.add_argument('-t', '--fit_type', default='e', help='data type: energy alone, energy & force, energy, force & hessian')
    args = parser.parse_args()
    extract_data(args.f_immm_input, args.atom_file,args.fit_type) 

if __name__ == '__main__':
    main()
