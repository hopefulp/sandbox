#!/home/joonho/anaconda3/bin/python

import argparse
import os
import sys


SPIN_LABELS = {
    0: 'up',
    1: 'down',
}


def _next_nonblank(lines, idx):
    while idx < len(lines) and not lines[idx].split():
        idx += 1
    return idx


def read_eigenval(fname):
    with open(fname, 'r') as f:
        lines = f.readlines()

    if len(lines) < 6:
        raise ValueError(f"{fname} is too short for EIGENVAL format")

    header = lines[5].split()
    if len(header) < 3:
        raise ValueError(f"cannot read nelect/nkpoints/nbands from line 6 of {fname}")

    nelect = int(float(header[0]))
    nkpoints = int(float(header[1]))
    nbands = int(float(header[2]))

    states = []
    idx = 6
    nspin = None

    for ikpt in range(1, nkpoints + 1):
        idx = _next_nonblank(lines, idx)
        if idx >= len(lines):
            raise ValueError(f"missing k-point {ikpt} block in {fname}")

        k_fields = lines[idx].split()
        if len(k_fields) < 4:
            raise ValueError(f"cannot read k-point {ikpt} line: {lines[idx].rstrip()}")
        kpoint = tuple(float(x) for x in k_fields[:3])
        weight = float(k_fields[3])
        idx += 1

        for _ in range(nbands):
            idx = _next_nonblank(lines, idx)
            if idx >= len(lines):
                raise ValueError(f"missing band rows after k-point {ikpt} in {fname}")

            fields = lines[idx].split()
            band = int(fields[0])

            if len(fields) == 3:
                row_nspin = 1
                energies = [float(fields[1])]
                occupations = [float(fields[2])]
            elif len(fields) >= 5:
                row_nspin = 2
                energies = [float(fields[1]), float(fields[2])]
                occupations = [float(fields[3]), float(fields[4])]
            else:
                raise ValueError(f"cannot read band row: {lines[idx].rstrip()}")

            if nspin is None:
                nspin = row_nspin
            elif nspin != row_nspin:
                raise ValueError(f"inconsistent spin columns near k-point {ikpt}, band {band}")

            for ispin, (energy, occupation) in enumerate(zip(energies, occupations)):
                states.append({
                    'energy': energy,
                    'occupation': occupation,
                    'spin': ispin,
                    'band': band,
                    'kpoint_index': ikpt,
                    'kpoint': kpoint,
                    'weight': weight,
                })

            idx += 1

    return {
        'nelect': nelect,
        'nkpoints': nkpoints,
        'nbands': nbands,
        'nspin': nspin or 1,
        'states': states,
    }


def find_vbm(eigen_data, occ_tol=1.0e-6, spin=None):
    occupied = []
    for state in eigen_data['states']:
        if spin is not None and state['spin'] != spin:
            continue
        if state['occupation'] > occ_tol:
            occupied.append(state)

    if not occupied:
        return None

    return max(occupied, key=lambda x: x['energy'])


def find_cbm(eigen_data, occ_tol=1.0e-6, spin=None):
    unoccupied = []
    for state in eigen_data['states']:
        if spin is not None and state['spin'] != spin:
            continue
        if state['occupation'] <= occ_tol:
            unoccupied.append(state)

    if not unoccupied:
        return None

    return min(unoccupied, key=lambda x: x['energy'])


def read_fermi_from_doscar(fname='DOSCAR'):
    if not os.path.exists(fname):
        return None

    with open(fname, 'r') as f:
        for _ in range(5):
            next(f, None)
        line = next(f, '').split()

    if len(line) < 4:
        return None
    return float(line[3])


def read_fermi_from_outcar(fname='OUTCAR'):
    if not os.path.exists(fname):
        return None

    efermi = None
    with open(fname, 'r') as f:
        for line in f:
            fields = line.split()
            if len(fields) >= 3 and fields[0] == 'E-fermi':
                efermi = float(fields[2])

    return efermi


def read_fermi(doscar='DOSCAR', outcar='OUTCAR'):
    efermi = read_fermi_from_doscar(doscar)
    if efermi is not None:
        return efermi, doscar

    efermi = read_fermi_from_outcar(outcar)
    if efermi is not None:
        return efermi, outcar

    return None, None


def print_band_edge(name, band_edge, label=None):
    if band_edge is None:
        if label:
            print(f"{name} ({label}): not found")
        else:
            print(f"{name}: not found")
        return

    prefix = name
    if label:
        prefix += f" ({label})"

    print(
        f"{prefix}: {band_edge['energy']:.8f} eV  "
        f"occ={band_edge['occupation']:.6f}  "
        f"kpt={band_edge['kpoint_index']} {band_edge['kpoint']}  "
        f"band={band_edge['band']}"
    )


def main():
    parser = argparse.ArgumentParser(description='Extract valence band maximum from VASP EIGENVAL')
    parser.add_argument('file', nargs='?', default='EIGENVAL', help='read EIGENVAL')
    parser.add_argument('-d', '--doscar', default='DOSCAR',
                        help='read Fermi level from DOSCAR')
    parser.add_argument('-f', '--fermi_file', default='OUTCAR',
                        help='fallback file to read Fermi level from OUTCAR')
    parser.add_argument('-o', '--occ_tol', type=float, default=1.0e-6,
                        help='occupation threshold for occupied states')
    args = parser.parse_args()

    try:
        eigen_data = read_eigenval(args.file)
    except (OSError, ValueError) as err:
        print(f"eigen.py: {err}", file=sys.stderr)
        sys.exit(1)

    print(
        f"NELECT={eigen_data['nelect']} "
        f"NKPTS={eigen_data['nkpoints']} "
        f"NBANDS={eigen_data['nbands']} "
        f"ISPIN={eigen_data['nspin']}"
    )
    efermi, fermi_source = read_fermi(args.doscar, args.fermi_file)
    if efermi is not None:
        print(f"E-fermi={efermi:.8f} eV ({fermi_source})")
    else:
        print(f"E-fermi: not found in {args.doscar} or {args.fermi_file}")

    print_band_edge('VBM', find_vbm(eigen_data, args.occ_tol))
    print_band_edge('CBM', find_cbm(eigen_data, args.occ_tol))
    if eigen_data['nspin'] == 2:
        for ispin in range(2):
            print_band_edge('VBM', find_vbm(eigen_data, args.occ_tol, spin=ispin), SPIN_LABELS[ispin])
            print_band_edge('CBM', find_cbm(eigen_data, args.occ_tol, spin=ispin), SPIN_LABELS[ispin])


if __name__ == '__main__':
    main()
