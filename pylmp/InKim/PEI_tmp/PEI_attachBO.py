#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys
import os
import random
import bgf
import bgftools
import LAMMPS_trj2bgf 
import nutils as nu
import tqdm

usage = "PEI_attachBO.py bgf_file bo_file out_file"
if len(sys.argv) < 2:
    print(usage)
    sys.exit(0)

curr_dir = os.path.abspath('.')

bgf_file = sys.argv[1]
bo_file = sys.argv[2]
out_file = sys.argv[3]

mybgf = bgf.BgfFile(bgf_file)

# find primary amines
pri_amines = []
for atom in mybgf.a:
    if atom.rName == "PRI" and "N" in atom.ffType and "N" in atom.aName:
        pri_amines.append(atom)

# attach
delatoms = []
for N_atom in tqdm.tqdm(pri_amines, ncols=80, desc="Analyzing PRI amines.. "):
    # chooise a random H to detach from the primary amine
    _ = [];
    for ano in N_atom.CONECT:
        atom = mybgf.getAtom(ano)
        if "H" in atom.ffType and "H" in atom.aName:
            _.append(ano)
    H_atom_body = mybgf.getAtom(random.choice(_))
    
    # buthyloxide
    new_bo = bgf.BgfFile(bo_file)

    # H15 will be removed and linked to the primary amine
    for atom in new_bo.a:
        if "H15" in atom.aName and "X" in atom.chain:
            H_atom_bo = atom
        elif "C" in atom.ffType and "C" in atom.aName and "H" in atom.chain:
            C_atom_bo = atom

    # move BO to link point
    (x, y, z) = (H_atom_body.x, H_atom_body.y, H_atom_body.z)
    for atom in new_bo.a:
        atom.x += x
        atom.y += y
        atom.z += z
        atom.rName = "BOX"
        atom.rNo = N_atom.rNo + 5000

    # merge
    mybgf = mybgf.merge(new_bo)
    mybgf.renumber()

    # connect
    N_atom.connect(C_atom_bo)
    N_atom.rName = "SEC"
    N_atom.charge += H_atom_body.charge
    C_atom_bo.charge += H_atom_bo.charge

    #delatoms = [mybgf.a2i[H_atom_body.aNo], mybgf.a2i[H_atom_bo.aNo]]
    delatoms.append(mybgf.a2i[H_atom_body.aNo])
    delatoms.append(mybgf.a2i[H_atom_bo.aNo])

mybgf.delAtoms(delatoms)
mybgf.renumber()

mybgf.CRYSTX = [50.0, 50.0, 50.0, 90.0, 90.0, 90.0]
mybgf.PERIOD = "111"
mybgf.AXES = "ZYX"
mybgf.SGNAME = "P 1                  1    1"
mybgf.CELLS = [-1, 1, -1, 1, -1, 1]

# align com
_ = bgftools.getCom(mybgf, '/home/noische/ff/DREIDING2.21.ff')
for atom in mybgf.a:
    atom.x = atom.x - _[0] + mybgf.CRYSTX[0]/2.0
    atom.y = atom.y - _[1] + mybgf.CRYSTX[1]/2.0
    atom.z = atom.z - _[2] + mybgf.CRYSTX[2]/2.0

# save
mybgf.saveBGF('temp.bgf')

# LAMMPS environments
mpi_command = "/opt/intel/Compiler/composer_xe_2011_sp1.13.367/mpi/openmpi-1.6.3/bin/mpirun -n 12 " # kdft
mpi_command = ""
lammps_command = "/qcfs/noische/program/kdft/lammps/bin/lmp_kdft "  # kdft

# minimization
nu.shutup()
createLammpsInput = "~tpascal/scripts/createLammpsInput.pl" + " -b temp.bgf " + " -f /home/noische/ff/DREIDING2.21.ff" + " -s temp " + " -o 'no shake' -t min " + " > /dev/null"
_ = os.system(createLammpsInput)
in_file = "in.temp"
data_file = "data.temp"
os.system("sed -i 's/dielectric      1/dielectric      72/' " + in_file)
os.system("sed -i 's/kspace_style    pppm 0.0001/kspace_style    none/' " + in_file)
os.system("sed -i 's/boundary        p p p/boundary        s s s/' " + in_file)
os.system("sed -i 's/1.0e-4 1.0e-4 500 5000/1.0e-6 1.0e-6 5000 50000/' " + in_file)
#os.system("sed -i 's/2.0 multi/5.0 multi/' " + in_file)
os.system("sed -i 's/lj\/charmm\/coul\/long\/opt 7.5 8.50000/lj\/cut\/coul\/debye 0.142 10/' " + in_file)
os.system("sed -i 's/every 2 delay 4/every 1 delay 0/' " + in_file)
#os.system("sed -i 's/0.000000  50.000000/-10.000000  10.000000/' " + data_file)
os.system("sed -i 's/0 # X/0 0 # X/' " + data_file)
os.system("sed -i 's/Impropers//' " + data_file)
runLammps = mpi_command + lammps_command + " -in in.temp -log lammps.log >> temp.log"
_ = os.system(runLammps)
LAMMPS_trj2bgf.getLAMMPSTrajectory("temp.bgf", "temp.min.lammpstrj", "temp.bgf", -1, False, True)

mybgf = bgf.BgfFile('temp.bgf')
nu.say()

mybgf.saveBGF(out_file)
