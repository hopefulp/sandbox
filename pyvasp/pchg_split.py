#!/home/joonho/anaconda3/bin/python

import argparse
import os
import sys


def is_int_line(line):
    fields = line.split()
    if not fields:
        return False
    try:
        [int(x) for x in fields]
    except ValueError:
        return False
    return True


def grid_size(line):
    fields = line.split()
    if len(fields) != 3:
        return None
    try:
        nx, ny, nz = [int(x) for x in fields]
    except ValueError:
        return None
    return nx * ny * nz


def get_header_and_grid(lines):
    if len(lines) < 8:
        raise ValueError("file is too short for a VASP PARCHG/CHGCAR format")

    count_idx = 5 if is_int_line(lines[5]) else 6
    if count_idx >= len(lines):
        raise ValueError("cannot find atom-count line")

    try:
        natom = sum(int(x) for x in lines[count_idx].split())
    except ValueError:
        raise ValueError(f"cannot read atom counts from line {count_idx + 1}")

    coord_idx = count_idx + 1
    if coord_idx < len(lines) and lines[coord_idx].strip().lower().startswith('s'):
        coord_idx += 1

    grid_idx = coord_idx + 1 + natom
    while grid_idx < len(lines) and not lines[grid_idx].split():
        grid_idx += 1

    if grid_idx >= len(lines):
        raise ValueError("cannot find charge-grid dimensions")

    ngrid = grid_size(lines[grid_idx])
    if ngrid is None:
        raise ValueError(f"cannot read charge-grid dimensions from line {grid_idx + 1}")

    return lines[:grid_idx + 1], grid_idx, ngrid


def parse_grid_value(field, line_no, field_no, allow_overflow, block_name):
    if set(field) == {'*'}:
        if allow_overflow:
            return None
        raise ValueError(
            f"VASP overflow value found in {block_name} grid "
            f"at line {line_no}, field {field_no}: {field}. "
            "The exact value is not recoverable from PARCHG; rerun VASP with an output "
            "that avoids overflow, or use --interpolate-overflow for an approximate split."
        )

    try:
        return float(field)
    except ValueError:
        raise ValueError(f"cannot read {block_name} grid value at line {line_no}, field {field_no}: {field}")


def interpolate_missing_values(values):
    missing = [i for i, value in enumerate(values) if value is None]
    if not missing:
        return values, 0

    fixed = values[:]
    nvalue = len(values)
    idx = 0
    while idx < nvalue:
        if fixed[idx] is not None:
            idx += 1
            continue

        start = idx
        while idx < nvalue and fixed[idx] is None:
            idx += 1
        end = idx - 1
        prev_idx = start - 1
        next_idx = idx

        if prev_idx >= 0 and next_idx < nvalue and fixed[prev_idx] is not None and fixed[next_idx] is not None:
            step = (fixed[next_idx] - fixed[prev_idx]) / (next_idx - prev_idx)
            for miss_idx in range(start, end + 1):
                fixed[miss_idx] = fixed[prev_idx] + step * (miss_idx - prev_idx)
        elif prev_idx >= 0 and fixed[prev_idx] is not None:
            for miss_idx in range(start, end + 1):
                fixed[miss_idx] = fixed[prev_idx]
        elif next_idx < nvalue and fixed[next_idx] is not None:
            for miss_idx in range(start, end + 1):
                fixed[miss_idx] = fixed[next_idx]
        else:
            raise ValueError("all charge-grid values are overflow markers; cannot interpolate")

    return fixed, len(missing)


def read_values(lines, start_idx, nvalue, allow_overflow=False, block_name='charge'):
    values = []
    idx = start_idx
    while idx < len(lines) and len(values) < nvalue:
        for field_no, field in enumerate(lines[idx].split(), start=1):
            values.append(parse_grid_value(field, idx + 1, field_no, allow_overflow, block_name))
            if len(values) == nvalue:
                break
        idx += 1

    if len(values) != nvalue:
        raise ValueError(f"expected {nvalue} grid values, found {len(values)}")

    values, ninterpolated = interpolate_missing_values(values)
    return values, idx, ninterpolated


def find_next_grid(lines, start_idx, ngrid):
    for idx in range(start_idx, len(lines)):
        if grid_size(lines[idx]) == ngrid:
            return idx
    return None


def read_ispin_from_outcar(fname='OUTCAR'):
    if not os.path.isfile(fname):
        return None

    ispin = None
    with open(fname, 'r') as f:
        for line in f:
            fields = line.split()
            if len(fields) >= 3 and fields[0] == 'ISPIN':
                try:
                    ispin = int(fields[2])
                except ValueError:
                    pass
    return ispin


def write_charge_file(fname, header, values, per_line=5):
    with open(fname, 'w') as f:
        f.writelines(header)
        for i in range(0, len(values), per_line):
            chunk = values[i:i + per_line]
            f.write(' '.join(f"{value:18.11E}" for value in chunk) + "\n")


def split_parchg(fname, up_file='PARCHG_up', down_file='PARCHG_down', allow_overflow=False):
    with open(fname, 'r') as f:
        lines = f.readlines()

    header, grid_idx, ngrid = get_header_and_grid(lines)
    total, next_idx, ntotal_interpolated = read_values(
        lines, grid_idx + 1, ngrid, allow_overflow, 'total charge'
    )
    mag_grid_idx = find_next_grid(lines, next_idx, ngrid)
    if mag_grid_idx is None:
        ispin = read_ispin_from_outcar()
        if ispin is not None:
            raise ValueError(
                f"spin/magnetization density block was not found; OUTCAR reports ISPIN={ispin}"
            )
        raise ValueError("spin/magnetization density block was not found; this PARCHG is not spin-polarized")

    mag, _, nmag_interpolated = read_values(
        lines, mag_grid_idx + 1, ngrid, allow_overflow, 'magnetization'
    )
    up = [(t + m) * 0.5 for t, m in zip(total, mag)]
    down = [(t - m) * 0.5 for t, m in zip(total, mag)]

    write_charge_file(up_file, header, up)
    write_charge_file(down_file, header, down)

    print(f"Grid points: {ngrid}")
    print(f"Spin-up charge: {up_file}")
    print(f"Spin-down charge: {down_file}")
    if ntotal_interpolated or nmag_interpolated:
        print(
            "Interpolated overflow values: "
            f"total={ntotal_interpolated}, magnetization={nmag_interpolated}"
        )


def main():
    parser = argparse.ArgumentParser(description='Split spin-polarized PARCHG into up and down charge densities')
    parser.add_argument('file', nargs='?', default='PARCHG', help='PARCHG file to read')
    parser.add_argument('-u', '--up', default='PARCHG_up', help='output spin-up PARCHG filename')
    parser.add_argument('-d', '--down', default='PARCHG_down', help='output spin-down PARCHG filename')
    parser.add_argument(
        '--interpolate-overflow',
        action='store_true',
        help='replace VASP *********** charge-grid overflow markers by linear interpolation',
    )
    args = parser.parse_args()

    try:
        split_parchg(args.file, args.up, args.down, args.interpolate_overflow)
    except (OSError, ValueError) as err:
        print(f"pchg_split.py: ERROR: {err}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
