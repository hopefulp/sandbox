#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os
import sys

key_match = {'GMAX': 'local part'}
nlines = {'local part': 1}


def _input_file(path, basename):
    return os.path.join(path, basename) if os.path.isdir(path) else path


def extract_zvals(potcar):
    """Read one ZVAL for each potential concatenated in POTCAR."""
    potcar = _input_file(potcar, 'POTCAR')
    zvals = []
    pattern = re.compile(r'\bZVAL\s*=\s*([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)')
    with open(potcar, 'r') as f:
        for line in f:
            match = pattern.search(line)
            if match:
                zvals.append(float(match.group(1)))
    if not zvals:
        raise ValueError(f"ZVAL was not found in {potcar}")
    return zvals


def read_poscar_composition(poscar):
    """Return species labels (when present) and atom counts from POSCAR."""
    poscar = _input_file(poscar, 'POSCAR')
    with open(poscar, 'r') as f:
        lines = f.readlines()
    if len(lines) < 7:
        raise ValueError(f"{poscar} is too short to be a POSCAR")

    fields = lines[5].split()
    try:
        counts = [int(value) for value in fields]
        species = [f"group{i + 1}" for i in range(len(counts))]  # VASP 4
    except ValueError:
        species = fields
        try:
            counts = [int(value) for value in lines[6].split()]  # VASP 5+
        except ValueError as exc:
            raise ValueError(f"cannot read atom counts from {poscar}") from exc

    if len(species) != len(counts):
        raise ValueError(
            f"species/count mismatch in {poscar}: "
            f"{len(species)} species and {len(counts)} counts"
        )
    return species, counts


def calculate_nelect(potcar, poscar):
    """Calculate neutral-system NELECT = sum(number of atoms * ZVAL)."""
    species, counts = read_poscar_composition(poscar)
    zvals = extract_zvals(potcar)
    if len(counts) != len(zvals):
        raise ValueError(
            "POTCAR/POSCAR species mismatch: "
            f"{len(zvals)} ZVAL values but {len(counts)} atom-count groups"
        )

    terms = [count * zval for count, zval in zip(counts, zvals)]
    for symbol, count, zval, subtotal in zip(species, counts, zvals, terms):
        print(f"{symbol:>4}: {count} * {zval:g} = {subtotal:g}")
    nelect = sum(terms)
    print(f"NELECT = {nelect:g}")
    return nelect

def fextract(fname, kw):
    if os.path.isdir(fname):
        fname += '/POTCAR'
    kwlines=[]
    kvs=[]
    kw = kw.upper()
    kw_orig = kw
    if kw == 'ZVAL':
        zvals = extract_zvals(fname)
        print(f"ZVAL: {[f'{value:g}' for value in zvals]}")
        return 0
    if kw in key_match.keys():
        kw = key_match[kw]
    ### READ POTCAR
    with open(fname, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):          # line has "\n"
        #print line,         # print writes its own "\n"
        if re.search(kw, line):
            if kw in nlines.keys():
                kwlines.append(lines[i+nlines[kw]].strip())
            else:
                kwlines.append(line.strip())
            print(line.strip())
    if not kwlines:
        print(f"there is no kw {kw} in {fname}")
        return 1
    else:
        for kwline in kwlines:
            lstr = re.split(r'[\s;]+', kwline)      #r'[;\s]+')
            #print(lstr)
            ### if key = value line
            if lstr[0] == kw:
                kvs.append(lstr[2])
            ### if only value line
            else:
                kvs.append(lstr[0])
    if kvs:
        print(f"{kvs}")
        if kw == 'ENMAX':
            enmax = list(map(float, kvs))
            print(f"{enmax}")
            default = round((max(enmax) * 1.3)/50)*50
            print(f"{kw}: {max(enmax)} roundoff with 30% for cell relaxation: {default}")
        else:
            if len(kvs) == 1:
                print(f"{kw_orig}: {kvs[0]}")
            else:
                print(f"{kw_orig}: {kvs}")
    return 0

def main():
    parser = argparse.ArgumentParser(description='READ POTCAR')
    parser.add_argument('file', nargs='?', default='POTCAR', help='read POTCAR')
    parser.add_argument('-p', '--poscar', help='read atom counts from POSCAR and calculate NELECT from POTCAR ZVAL')
    parser.add_argument('-k', '--kw', default='enmax', help='extract value using key such as enmax,zval, gmax')
    parser.add_argument('-u', '--usage', action='store_true', help='print usage')
    args = parser.parse_args()

    if args.usage:
        print(f"potcar_kw.py POTCAR_H -k gmax 'only for local part'\
            \n\tdata_treat.py -i POTCAR2.H/0/POTCAR POTCAR2.H/0.1/POTCAR -l 57 25\
            \npotcar_kw.py POTCAR_H -k enmax\
            \npotcar_kw.py POTCAR -p POSCAR  # calculate NELECT from ZVAL and atom counts\
            ")
        sys.exit(0)
    
    try:
        if args.poscar:
            calculate_nelect(args.file, args.poscar)
        else:
            fextract(args.file, args.kw)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))

if __name__ == '__main__':
    main()
