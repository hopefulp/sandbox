#!/home/joonho/anaconda3/bin/python


import os
import numpy as np

from ase import Atoms, Atom
from ase.calculators.vasp import Vasp
from ase.lattice import bulk

from ase.optimize import QuasiNewton
from ase.visualize import view
from ase.io import write, read

from ase.structure import molecule
from ase.constraints import FixAtoms


def __main__():
    bulk = read('POSCAR') 

    #c = FixAtoms(mask=[s == 'Cu' for s in adsorbedH.get_chemical_symbols()])
    #adsorbedH.set_constraint(c)
    calc = Vasp(
            encut = 500,
    #        ispin = 2,
            prec = "Accurate",
    #		ismear = -1,
            sigma = 0.1,
            ediff = 1e-8,
            ediffg = -0.01,
            algo = "Fast",
            gga = "RP",
            xc = "PBE",
            kpts = (1,1,1),	
    #		isif = 0,
    #		ibrion = 5,
    #		nsw = 0,
    #		nfree = 2
    )

            
    adsorbedH.set_calculator(calc)
    electronicenergy = adsorbedH.get_potential_energy()

    print("electronic energy is %.5f"%electronicenergy)

    vib = Vibrations(adsorbedH, indices=[23], delta=0.01, nfree=2)
    vib.run()
    print(vib.get_frequencies())
    vib.summary()
    print(vib.get_mode(-1))
    vib.write_mode(-1)
    vib_energies = vib.get_energies()

    thermo = HarmonicThermo(vib_energies=vib_energies, electronicenergy=electronicenergy)
    thermo.get_entropy(temperature=298.15)
    thermo.get_internal_energy(temperature=298.15)
    thermo.get_free_energy(temperature=298.15)


    #exitcode =  os.system('vasp')
    exitcode = os.system('mpirun -np 4 vasp')

   
if __name__ == 'main':
    main()
