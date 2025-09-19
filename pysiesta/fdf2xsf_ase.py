from ase import Atoms
from ase.io import write
import numpy as np

def read_fdf_atoms(filename):
    species_map = {}   # species_index -> atomic_number
    symbols = []
    positions = []

    with open(filename) as f:
        lines = f.readlines()

    # Parse ChemicalSpeciesLabel
    in_block = False
    for line in lines:
        if line.strip().lower().startswith("%block chemicalspecieslabel"):
            in_block = True
            continue
        if line.strip().lower().startswith("%endblock chemicalspecieslabel"):
            in_block = False
            continue
        if in_block:
            parts = line.split()
            if len(parts) >= 2:
                idx = int(parts[0])
                Z = int(parts[1])
                species_map[idx] = Z

    # Parse AtomicCoordinatesAndAtomicSpecies
    in_block = False
    for line in lines:
        if line.strip().lower().startswith("%block atomiccoordinatesandatomicspecies"):
            in_block = True
            continue
        if line.strip().lower().startswith("%endblock atomiccoordinatesandatomicspecies"):
            in_block = False
            continue
        if in_block:
            parts = line.split()
            if len(parts) >= 4:
                x, y, z = map(float, parts[:3])
                idx = int(parts[3])
                Z = species_map[idx]
                symbols.append(Z)
                positions.append([x, y, z])

    return Atoms(numbers=symbols, positions=np.array(positions))

# --- Usage ---
atoms = read_fdf_atoms("h2o.fdf")
atoms.center(vacuum=3.0)  # optional: add 3 Å vacuum around
atoms.write("h2o.xsf")    # save for XCrySDen / ASE
