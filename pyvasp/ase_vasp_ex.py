#!/home/joonho/anaconda3/bin/python
from ase.build import bulk
from ase.calculators.vasp import Vasp

si = bulk('Si')

mydir = 'bandstructure'    # Directory where we will do the calculations

# Make self-consistent ground state
calc = Vasp(kpts=(4, 4, 4), directory=mydir)

si.calc = calc
si.get_potential_energy()  # Run the calculation

# Non-SC calculation along band path
kpts = {'path': 'WGX',     # The BS path
        'npoints': 30}     # Number of points along the path

calc.set(isym=0,           # Turn off kpoint symmetry reduction
         icharg=11,        # Non-SC calculation
         kpts=kpts)

# Run the calculation
si.get_potential_energy()

