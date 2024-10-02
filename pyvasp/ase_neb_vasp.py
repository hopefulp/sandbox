#!/home/joonho/anaconda3/bin/python

'''
for ASE-NEB
'''
import argparse
import os
import numpy as np
import sys, shutil, copy

### ASE
import ase
from ase.io import read
from ase.constraints import FixAtoms
#from ase.build.structure import molecule
from ase.neb import SingleCalculatorNEB
from ase.optimize import BFGS   # QuasiNewton
from ase.neb import NEBTools

import matplotlib.pyplot as plt

### VASP
from ase import Atoms, Atom
from ase.calculators.vasp import Vasp
from ase.units import eV
#from ase.calculators.siesta.parameters import Species, PAOBasisBlock
'''
env
    export VASP_PP_PATH=/home/joonho/sandbox/pyvasp/POTCARfe
    export VASP_COMMAND=mpirun -np 4 /TGM/Apps/VASP/OLD_BIN/5.4.4/O2/NORMAL/vasp_bin/vasp.5.4.4.pl2.O2.NORMAL.std.x
'''


def run_neb(ini_geo, fin_geo, n_images, npar):
    '''
    ini_geo initial atomic structure
    fin_geo final   atomic structure
    nimages without 00 and 0f
    npar    parallel for wavefunction
    '''

    ### Step 1 : Set SIESTA calculator
    print(ase.__version__)
    ase_calculator  = Vasp(
        label='neb-vasp',
        encut = 500,
        prec = 'Accurate',
        algo = 'Fast',
        sigma = 0.05,   # if metal atoms, SIGMA= 0.1 / 0.2(default) (ex. Transition metals), if semi/insulator (SIGMA= 0.05)
        ediff = 1e-4,
        ediffg = -0.01,
        gga='PE',
        kpts = (2, 4, 1)
        )
    ### Step 2 : Structures
    struct_ini = read(ini_geo) 
    struct_fin = read(fin_geo) 

    ### Step 3 : Optimize initial and fianl end-points

    struct_ini.set_calculator(copy.deepcopy(ase_calculator))
    qn = BFGS(struct_ini, trajectory='initial.traj', logfile='init.log')
    qn.run(fmax=0.05)
    struct_fin.set_calculator(copy.deepcopy(ase_calculator))
    qn = BFGS(struct_fin, trajectory='final.traj', logfile ='fin.log')
    qn.run(fmax=0.05)

    ### Step 4 : NEB using ASE

    initial_ase = read('initial.traj')
    final_ase = read('final.traj')
    constraint = FixAtoms(mask=[atom.tag > 1 for atom in initial_ase])

    images_ase = [initial_ase]
    for i in range(1, n_images-1):
        image = initial_ase.copy()
        image.set_calculator(ase_calculator)
        image.set_constraint(constraint)
        images_ase.append(image)

    images_ase.append(final_ase)

    neb_ase = SingleCalculatorNEB(images_ase, climb=True)
    neb_ase.interpolate()

    qn_ase = BFGS(neb_ase, trajectory='neb_ase.traj')
    qn_ase.run(fmax=0.05)

    ### Step 5 : Summary of results

    # NEB ASE
    print('\nSummary of the results: \n')

    atoms_ase = read('neb_ase.traj', ':')
    n_eval_ase = int(len(atoms_ase) - 2 * (len(atoms_ase)/n_images))

    print('Number of function evaluations CI-NEB implemented in ASE:', n_eval_ase)

    # Plot ASE NEB:
    nebtools_ase = NEBTools(images_ase)

    Sf_ase = nebtools_ase.get_fit()[2]
    Ef_ase = nebtools_ase.get_fit()[3]

    Ef_neb_ase, dE_neb_ase = nebtools_ase.get_barrier(fit=False)
    nebtools_ase.plot_band()

    plt.savefig('neb_result.png', dpi=300, bbox_inches ="tight")

def main():
    parser = argparse.ArgumentParser(description='run ASE NEB using only POSCARs')
    parser.add_argument('struct_ini',  help='initial geometry for siesta')
    parser.add_argument('struct_fin',  help='final geometry for siesta')
    parser.add_argument('-n', '--nimages', default=7, help='number of calculation images except 00 and last one')
    parser.add_argument('-p', '--npar', default=4, help='npar for vasp parallel')
    parser.add_argument('-u', '--usage', action='store_true', help='Usage')
    args = parser.parse_args()	 

    if args.usage:
        print(f"Set env::\
                \n\texport VASP_PP_PATH=/home/joonho/sandbox/pyvasp/POTCARfe\
                \n\texport VASP_COMMAND=mpirun -np 4 /TGM/Apps/VASP/OLD_BIN/5.4.4/O2/NORMAL/vasp_bin/vasp.5.4.4.pl2.O2.NORMAL.std.x\
                Run::\
                \n\tCLI-python\
                \n\t    ase_neb_vasp.py POSCAR.HfSe2mO2sglini POSCAR.HfSe2mO2sglfin -n 7\
                \n\tPBS\
                \n\t    qsub ase_vasp.sh\
                ")
        sys.exit(0)

    run_neb(args.struct_ini, args.struct_fin, args.nimages, args.npar)

if __name__ == '__main__':
    main()

