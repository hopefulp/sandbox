from ase.calculators.vasp import Vasp
from ase.lattice.surface import *
from ase.visualize import view
from ase import Atoms
from ase.io import write
from ase.io import read 
from ase.thermochemistry import HarmonicThermo 
from ase.vibrations import Vibrations
from ase.structure import molecule
from ase.constraints import FixAtoms


adsorbedH = read('CONTCAR') 

#c = FixAtoms(mask=[s == 'Cu' for s in adsorbedH.get_chemical_symbols()])
#adsorbedH.set_constraint(c)

calc = Vasp(
		encut = 500,
		ispin = 2,
		prec = "Accurate",
#		ismear = -1,
		sigma = 0.1,
		ediff = 1e-8,
		ediffg = -0.01,
		algo = "Fast",
		gga = "RP",
		xc = "PBE",
		gamma = 1,
		#kpts = (2,2,2),	
#		isif = 0,
#		ibrion = 5,
#		nsw = 0,
#		nfree = 2
)

		
adsorbedH.set_calculator(calc)
electronicenergy = adsorbedH.get_potential_energy()

print "electronic energy is %.5f"%electronicenergy

vib = Vibrations(adsorbedH, indices=[13,14], delta=0.01, nfree=2)
vib.run()
print vib.get_frequencies()
vib.summary()
print vib.get_mode(-1)
vib.write_mode(-1)
vib_energies = vib.get_energies()

thermo = HarmonicThermo(vib_energies=vib_energies, electronicenergy=electronicenergy)
thermo.get_entropy(temperature=298.15)
thermo.get_free_energy(temperature=298.15)
