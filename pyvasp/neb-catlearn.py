#!/home/joonho/anaconda3/bin/python
import sys, shutil, copy
from ase.io import read
from ase.optimize import BFGS
from ase.calculators.vasp import Vasp
from catlearn.optimize.mlneb import MLNEB
from ase.neb import NEBTools
from catlearn.optimize.tools import plotneb

### Read input files
struct_init = read('CONTCAR_initial')
struct_fin  = read('CONTCAR_final')

### Set calculator
ase_calculator = Vasp(encut=400,
                      xc='PBE',
                      gga='PE',
                      istart = 0,
                      lwave=False,
                      lcharg=False,
                      kpts = (3, 3, 1),
                      ediffg=-0.01,
                      ediff=1e-6,
                      ibrion=1,
                      nsw=20,
                      ismear=1,
                      sigma=0.20,
                      algo='Normal',
                      prec='Normal'
                      )

# Optimize initial state:
struct_init.set_calculator(copy.deepcopy(ase_calculator))
qn = BFGS(struct_init, trajectory='initial.traj')
qn.run(fmax=0.01)

# Optimize final state:
struct_fin.set_calculator(copy.deepcopy(ase_calculator))
qn = BFGS(struct_fin, trajectory='final.traj')
qn.run(fmax=0.01)

# CatLearn NEB:

neb_catlearn = MLNEB(start='initial.traj', end='final.traj',
                    ase_calc=copy.deepcopy(ase_calculator),
                    n_images=11,
                    )

neb_catlearn.run(fmax=0.05, trajectory='ML-NEB.traj')

