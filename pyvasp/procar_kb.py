#!/home/joonho/anaconda3/envs/pmg39/bin/python
'''
conversion of procar_kb.pl to procar_kb.py 2025.08.04 by Joonho Park
'''

import argparse
import os
import sys
from libdoscar import l_map
import numpy as np

def get_energy_range(energy_str):
    if len(energy_str) == 2:
        Emin = energy_str[0]
        Emax = energy_str[1]
    else:
        E_pivot = float(energy_str[0])
        delta = 0.5
        Emin = E_pivot - delta
        Emax = E_pivot + delta
    return Emin, Emax

'''
def orbital_columns(orbital):
    if orbital == "s":
        return [1]
    elif orbital == "p":
        return list(range(2, 5))  # px, py, pz
    elif orbital == "d":
        return list(range(5, 10))  # dxy, dyz, dz2, dxz, dx2
    return [1, 2, 3, 4, 5, 6, 7, 8, 9]  # fallback: all orbitals
'''

def sum_orbitals(line: str, col_indices: list[int]) -> float:
    values = list(map(float, line.split()))
    #orbital_values = values[1:-1]  # exclude atom index and total
    return sum(values[i] for i in col_indices)          # use col 1(s) ~ 10 including total



def analyze_procar(fin, Emin, Emax, atom_indices, orbitals, density_crit, poptions):
    '''
    In  fin             PROCAR
        Emin, Emax      Energy range
        atom_indices    atom indices to check which starts from 1
        orbitals        each atom's orbital list
        density_crit    control parameter to check density
        doptions        display options
            a           all if k,b above density_crit
            t:int           display list for top 10
    '''
    with open(fin, 'r') as f:
        lines = f.readlines()

    print(f"Reading: {fin}")

    kpoints = bands = atoms = None
    energy = None
    k_b_dos = []
    collecting = False

    # use different l's for each atom
    orbital_cols2D = []
    for i, atom in enumerate(atom_indices):
        if orbitals[i] == 't':
            orbitals_cols = [10]
        else:
            orbitals_cols - l_map[orbitals[i]]
        orbital_cols2D.append(orbitals_cols)
    #orbital_cols2D = orbital_columns(orbitals)
    max_atom_block = None

    ### PROCAR format
    # 0-th: procar lm decomposed
    # 1-st: kpoints, nbands, nions
    # 2-nd: blank                   3 lines
    # kpoints: 2 lines
    # band: 2 lines
    # l lead: 1 lines
    # end of ion block: total+blank: 2 lines

    parts = lines[1].split()
    nkpoints = int(parts[3])
    nbands = int(parts[7])
    natoms = int(parts[11])
    weight = 1.0 / nkpoints
    #max_atom_block = atoms + 4
    print(f"nkpoints {nkpoints} nbands {nbands} natoms {natoms}")

    i = 3
    while i < len(lines):
        ### inside while-loop, i should be treated explicitely
        ### k-point with 2 lines
        #print(f"i {i}: line {lines[i]}")
        if 'k-point' in lines[i] or 'k-point' in lines[i+1]:    # in case new k-point, 1 more blank line
            if 'k-point' in lines[i]:
                pass
            elif 'k-point' in lines[i+1]:
                i += 1
            parts = lines[i].strip().split()
            ik_point = parts[1]
            i += 2                              # k-points line has 1 blank line
            continue
        ### band with 2 lines : in k, in band if energy in range, go on
        if "band" in lines[i]:
            parts = lines[i].split()
            iband  = int(parts[1])
            energy = float(parts[4])
            collecting = Emin <= energy <= Emax          # collecting gets True or False in E-range
            skip_counter = 0
            i += 2                              # band has 1 blank line & ion head goes into ion-block
            continue
        if collecting:
            dos_vals=[]
            ion_block = lines[i+1:i+1+natoms]           # ion_block without head
            for j, iatom in enumerate(atom_indices):    # j: ordering in input atom_indices
                parts = ion_block[iatom-1].split()
                if iatom != int(parts[0]):
                    print(f"error in format in ion block: input {iatom} line {parts[0]} {i}-step")
                    sys.exit(10)
                #print(f"ion-line {ion_block[iatom-1]}, orbital columns {orbital_cols2D[j]}")
                dos_val = sum_orbitals(ion_block[iatom-1], orbital_cols2D[j])
                if dos_val >= density_crit:
                    dos_vals.append(dos_val)
            if len(dos_vals) == len(atom_indices):
                #k_b_dos.append((energy, ik_point, iband), dos_vals))        # if all atoms over criterion
                lenekb = [energy, ik_point, iband]
                lenekb.extend(dos_vals)        # if all atoms over criterion
                k_b_dos.append(lenekb)
        i += natoms + 3                                 # ion head and total/blank - all 3 lines
        #if i > 1000:
        #    sys.exit(100)
        #print(f"{i}-th line in PROCAR loop")

    if not k_b_dos:
        print("No significant DOS found in energy window.")
        return

    if 't' in poptions:
        nsort = int(poptions[1:])
        for k_b_1dos in k_b_dos:
            #print(f"dos {k_b_1dos[-2]} {k_b_1dos[-1]}")
            dos_prod = k_b_1dos[-2] * k_b_1dos[-1]
            k_b_1dos.append(dos_prod)
        sorted_dos = sorted(k_b_dos, key=lambda row: row[-1], reverse=True)

        for row in sorted_dos:
            del row[-1]
        k_b_dos = sorted_dos[:nsort]
        

    ### print
    print(f"\nEnergy\t\tKpoint\tband", end='')
    for iatom in atom_indices:
        print(f"\tDOS[{iatom}]", end='')
    print("")
    
    total_dos = [0, 0]
    for e, ikpoint, iband, dos1, dos2 in k_b_dos:
        print(f"{e:.5f}\t{ikpoint}\t{iband}\t{dos1:.3f}\t{dos2:.3f}")
        total_dos[0] += dos1
        total_dos[1] += dos2

    print("\nSummed DOS (weighted):\t", end='     ')
    for i, atom in enumerate(atom_indices):
        print(f"{total_dos[i] * weight:8.3f}", end='')
    print("\n")


def main():
    parser = argparse.ArgumentParser(description="Analyze PROCAR for orbital overlap or DOS in a given energy range")
    parser.add_argument("directory", default=".", help="Directory containing PROCAR file")
    parser.add_argument("-e", "--energy", nargs='+', type=float, help="Energy pivot or range as Emin:Emax")
    parser.add_argument("-ia", "--atom_indices", nargs="+", type=int, help="List of atom indices (1-based)")
    parser.add_argument("-l", "--orbitals", nargs='*', default=['t','t'], choices=["s", "p", "d"], help="Angular momentum (orbital) to analyze")
    parser.add_argument("-c", "--density_crit", default=0.005, help="density criteria to check")
    parser.add_argument("-o", "--options", default='a', help="display option: 'a'll, 't'op 10")

    args = parser.parse_args()

    Emin, Emax = get_energy_range(args.energy)
    procar_path = os.path.join(args.directory, "PROCAR")

    if not os.path.isfile(procar_path):
        print(f"Error: PROCAR not found at {procar_path}")
        exit(1)
    else:
        print(f"read PROCAR {procar_path}")

    ### use the same number of l's as number of atoms
    if len(args.orbitals) < len(args.atom_indices):
        while len(args.orbitals) < len(args.atom_indices):
            args.ortibals.append('t')
    

    analyze_procar(procar_path, Emin, Emax, args.atom_indices, args.orbitals, args.density_crit, args.options)


if __name__ == "__main__":
    main()
