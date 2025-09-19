#!/home/joonho/anaconda3/bin/python

#!/usr/bin/env python3
import argparse
import numpy as np
from ase import Atoms
from ase.data import atomic_numbers
from libsci import generate_distances
import os, sys

def generate_diatomic(atom1, atom2, r):
    """Create ASE Atoms object for diatomic molecule with distance r (Å)."""
    atoms = Atoms(f'{atom1}{atom2}',
                  positions=[[0, 0, 0], [0, 0, r]],
                  pbc=False)
    atoms.center(vacuum=10.0)  # add vacuum padding
    return atoms

def get_num_species(atoms):
    atom_kinds = []
    for atom in atoms:
        if atom.symbol not in atom_kinds:
            atom_kinds.append(atom.symbol)
    print(f"Unique atoms species: {atom_kinds}")
    return len(atom_kinds)


def write_fdf(filename, atoms):
    """
    Write minimal SIESTA .fdf input structure from an ASE Atoms object.
    Includes    NumberOfAtoms
                NumberOfSpecies
                LatticeVectors
                ChemicalSpeciesLabel
                AtomicCoordinatesAndAtomicSpecies.
    """

    # --- Unique species dictionary (sorted for reproducibility) ---
    symbols = sorted(set([atom.symbol for atom in atoms]))
    species_dict = {sym: i+1 for i, sym in enumerate(symbols)}  # 1-based index

    with open(filename, "w") as f:
        # Header
        f.write(f"NumberOfAtoms {len(atoms)}\n")
        f.write(f"NumberOfSpecies {len(species_dict)}\n\n")


        # ChemicalSpeciesLabel block
        f.write("%block ChemicalSpeciesLabel\n")
        for sym, idx in species_dict.items():
            Z = atomic_numbers[sym]  # atomic number
            f.write(f"  {idx}  {Z}  {sym}\n")
        f.write("%endblock ChemicalSpeciesLabel\n\n")


        f.write("LatticeConstant 1.00000 Ang\n")

        # Lattice vectors (angstrom)
        f.write("%block LatticeVectors\n")
        for vec in atoms.cell:  # 3x3 matrix
            f.write(f"  {vec[0]:12.6f} {vec[1]:12.6f} {vec[2]:12.6f}\n")
        f.write("%endblock LatticeVectors\n\n")

        f.write("AtomicCoordinatesFormat Ang\n")

        # Atomic coordinates (Cartesian, Angstrom)
        f.write("%block AtomicCoordinatesAndAtomicSpecies\n")
        for atom in atoms:
            x, y, z = atom.position
            sym = atom.symbol
            idx = species_dict[sym]
            f.write(f"  {x:12.6f} {y:12.6f} {z:12.6f}  {idx}\n")
        f.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")

def main():
    parser = argparse.ArgumentParser(
        description="Generate bond dissociation geometries for diatomic molecule."
    )
    parser.add_argument("-a","--atoms", nargs=2, default=["H", "H"], help="Two atoms for dissociation)")
    parser.add_argument("-rmme","--rminmaxeq", nargs='*', type=float, default=[0.5, 10, 0.74], help="minimum maximum equilibrium distance in Å (default=0.5)")
    parser.add_argument("-np", "--npoints", type=int, default=10,  help="Number of geometries (default=10) between r_eq ~ r_max")
    parser.add_argument("-npl", "--npoints_left", type=int, default=3,  help="Number of geometries between r_min ~ r_eq")
    parser.add_argument("-o","--outdir", default="geometries", help="Output directory (default=geometries)")
    parser.add_argument("-p", "--prefix", default="H2p", help="File prefix for outputs (default=mol)")
    parser.add_argument('-u', '--usage',   action='store_true', help='print usage')

    args = parser.parse_args()

    if args.usage:
        print('Make two atom dissociation\
              \n\tbuild_struct.py -a H H -rmme 0.5 10 0.74 -np 15 -p H2 for neutral\
              \n\tbuild_struct.py -a H H -rmme 0.7 10 1.06 -np 15 -p H2p for H2 cation\
              ')
        sys.exit(0)

    os.makedirs(args.outdir, exist_ok=True)

    # Generate evenly spaced bond lengths
    atom1 = args.atoms[0]
    if len(args.atoms) == 2:
        atom2 = args.atoms[1]
    else:
        atom2 = atom1

    rmin = args.rminmaxeq[0]
    if len(args.rminmaxeq) >= 2:
        rmax = args.rminmaxeq[1]
    else:
        rmax = rmin + 10.0
    if len(args.rminmaxeq) == 3:
        req  = args.rminmaxeq[2]
    else:
        req = None

    if req:
        distances = generate_distances(rmin, req, rmax, n_short=args.npoints_left, n_long=args.npoints)
    else:
        distances = np.linspace(rmin, rmax, args.npoints)

    for i, r in enumerate(distances, 1):
        atoms = generate_diatomic(atom1, atom2, r)
        fname = os.path.join(args.outdir, f"{args.prefix}_{r:05.2f}.fdf")
        write_fdf(fname, atoms)
        print(f"Written: {fname} at r = {r:.3f} Å")

if __name__ == "__main__":
    main()
