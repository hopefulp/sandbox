#!/usr/bin/env python3
import os, re, sys
import argparse
import pandas as pd
from mplot2D import mplot_nvector

# Keywords to extract
ENERGY_KEYS = [
    "Kinetic",
    "Ion-electron",
    "Hartree",
    "Exch.-corr.",
    "Ion-ion",
    "Total"
]

def parse_stdout(file_path):
    """Extract energy values from a siesta stdout file."""
    values = {}
    try:
        with open(file_path, "r") as f:
            for line in f:
                for key in ENERGY_KEYS:
                    if key in line:
                        try:
                            values[key] = float(line.split()[-1])
                        except ValueError:
                            values[key] = None
    except FileNotFoundError:
        return None
    return values if values else None

def extract_distance(subdir):
    """Extract numeric distance from directory name."""
    match = re.search(r"([0-9]+\.?[0-9]*)", subdir)
    if match:
        return float(match.group(1))
    else:
        return subdir  # fallback: keep string


def collect_data(rootdir, wdir, siesta_outf):
    """ Search subdirs for stdout.txt and collect energies.
        if wdir: DFT runs under wdir under subdir: subdir/wdir
    """
    data = {}
    for subdir in os.listdir(rootdir):
        dft_dir = os.path.join(rootdir, subdir, wdir) if wdir else os.path.join(rootdir, subdir)
        if os.path.isdir(dft_dir) and siesta_outf in os.listdir(dft_dir):
            fpath = os.path.join(dft_dir, siesta_outf)
            ### make column label
            label = extract_distance(subdir)  # use subdir name as column label
            values = parse_stdout(fpath)
            if values:
                data[label] = values
    return data


def main():
    parser = argparse.ArgumentParser(
        description="Collect SIESTA energies from stdout.txt in subdirs"
    )
    parser.add_argument("rootdir", help="Root directory to search")
    parser.add_argument("-wd","--work_dir", default='1.dft', help="if there is a specific work_dir under sub_dir")
    parser.add_argument("-so","--siesta_output", default='stdout.txt', help="Siesta output filename")
    parser.add_argument(
        "-o", "--output", default="energies.csv",
        help="Output CSV file (default: energies.csv)"
    )
    args = parser.parse_args()

    data = collect_data(args.rootdir, args.work_dir, args.siesta_output)
    if not data:
        print("No stdout.txt files found.")
        return 1

    # Build DataFrame: rows = ENERGY_KEYS, columns = subdir names
    df = pd.DataFrame(data).reindex(ENERGY_KEYS)
    print(df)

    # Save CSV
    df.to_csv(args.output, float_format="%.6f")
    print(f"\nSaved table to {args.output}")

    ### how to extract x and y for plotting
    x = df.columns.to_numpy(dtype=float)              # (m,)

    ### selected rows
    selected_rows = ["Kinetic","Total"]
    y_total = df.loc["Total"].to_numpy(dtype=float)   # (m,)
    y_all = df.loc[selected_rows].to_numpy()    # (n, m)

    plot_dict = {}
    plot_dict['legends'] = selected_rows
    mplot_nvector(x, y_all, plot_dict=plot_dict)

    return 0



if __name__ == "__main__":
    main()
