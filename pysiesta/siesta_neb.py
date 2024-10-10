#!/usr/bin/python3

### Update : 02/04/21  Min Jong Noh & Hyeonwoo Yeo
###          09/30/24
'''
for ASE-NEB
'''
import argparse
import sys, shutil, copy
### ASE
import ase
from ase.build import fcc100, add_adsorbate
from ase.io import read
from ase.constraints import FixAtoms
from ase.neb import SingleCalculatorNEB
from ase.optimize import BFGS
import matplotlib.pyplot as plt
from ase.neb import NEBTools
### SIESTA
from ase import Atoms
from ase.calculators.siesta import Siesta
from ase.units import Ry, eV
from ase.calculators.siesta.parameters import Species, PAOBasisBlock

### Step 1 : Set SIESTA calculator
print(ase.__version__)
ase_calculator  = Siesta(\
    label='neb-siesta',
    xc='PBE',
    pseudo_qualifier='gga',
    mesh_cutoff=200 * Ry,
    energy_shift=160 * 0.001 * eV,
    basis_set='DZP',
    kpts=[6, 6, 1],
    fdf_arguments={'DM.MixingWeight': 0.05,
        'MaxSCFIterations': 700,
        'SCF.DM.Tolerance': 1e-4,
        'SCF.DM.Converge': True,
        'PAO.BasisType': 'split',
        'SlabDipoleCorrection': '.true.' ,
        'DM.NumberPulay': '10',
        'ElectronicTemperature': '300.0 K',
        'LatticeConstant': '1.0 Ang',
        '%block' :'LatticeVectors',
        '10.668566949    -6.159500000':'    -0.000000000',
        '0.000000000    12.319000000':'    -0.000000000',
        '0.000000000     0.000000000':'    46.179000000',
        '%endblock':' LatticeVectors'})

    def run_neb(ini_geo, fin_geo, n_images):
        '''
        ini_geo initial atomic structure
        fin_geo final   atomic structure
        nimages without 00 and 0f
        '''

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
        parser = argparse.ArgumentParser(description='prepare vasp input files: -s for POSCAR -p POTCAR -k KPOINTS and -i INCAR')
        parser.add_argument('--struct_ini',  help='initial geometry for siesta')
        parser.add_argument('--struct_fin',  help='final geometry for siesta')
        parser.add_argument('-n', '--nimages', default=7, help='number of calculation images except 00 and last one')
        args = parser.parse_args()	 

        run_neb(args.struct_ini, args.struct_fin, arga.nimages)

    if __name__ == '__main__':
        main()


