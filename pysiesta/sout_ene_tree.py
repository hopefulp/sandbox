#!/usr/bin/env python3
import os, re, sys
import argparse
import pandas as pd
from mplot2D import mplot_nvector

'''
to use PANDAS use dict of dicts
e.g. column-wise indexing
    data = {
        "H2p_00.00": {
            "Kinetic": 14.719166,
            "Ion-electron": -26.381439,
            "Hartree": 4.020258,
            "Exch.-corr.": -6.666122,
            "Ion-ion": -1.884185,
            "Total": -16.192322,
            "distance": 0.00
        },
        "H2p_00.70": {
            "Kinetic": 14.650230,
            "Ion-electron": -26.450002,
            "Hartree": 4.018332,
            ... 
        }
        ... 
    }
'''
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
        fpath = os.path.join(dft_dir, siesta_outf)
        if os.path.isfile(fpath):
            ### make column label
            #label = os.path.basename(subdir)[-5:]  # use subdir name as column label
            values = parse_stdout(fpath)
            if values:
                ### column label is subdir name
                label = os.path.basename(subdir)
                data[label] = values
                ### add extra distance rwo (last 5 chars parse as float)
                try:
                    dist_val = float(label[-5:])
                except ValueError:
                    dist_val = None
                data[label]["distance"] = dist_val
    return data


def main():
    parser = argparse.ArgumentParser(
        description="Collect SIESTA energies from stdout.txt in subdirs"
    )
    parser.add_argument("rootdir", nargs='?', default='cal', help="Root directory to search")
    parser.add_argument("-wd","--work_dir", default='1.dft', help="if there is a specific work_dir under sub_dir")
    parser.add_argument("-so","--siesta_output", default='stdout.txt', help="Siesta output filename")
    parser.add_argument("-sr","--selected_rows", nargs='+', default=['Ion-electron'], \
        choices=["Kinetic", "Ion-electron", "Hartree", "Exch.-corr.", "Ion-ion", "Total",'0','1','2','3','4','5'], help="row selection from energy table")
    parser.add_argument(
        "-o", "--output", default="energies.csv",
        help="Output CSV file (default: energies.csv)"
    )
    parser.add_argument('-u', '--usage',   action='store_true', help='print usage')

    args = parser.parse_args()

    if args.usage:
        print('Prepare Siesta energy table and plot by rows\
              \n\t$ sout_ene.py cal\
              \n\tFor plot:\
              \n\t    -sr --selectred_rows in ["Kinetic", "Ion-electron", "Hartree", "Exch.-corr.", "Ion-ion", "Total"] or index\
              \n\t$ sout_ene.py cal -sr 4\
              ')
        sys.exit(0)

    data = collect_data(args.rootdir, args.work_dir, args.siesta_output)
    if not data:
        print("No stdout.txt files found.")
        return 1

    ### add "distance" to ENERGY_KEYS for row order
    row_order = ENERGY_KEYS + ["distance"]
    # Build DataFrame: rows = ENERGY_KEYS, columns = subdir names
    df = pd.DataFrame(data).reindex(row_order)
    print(df)

    # Save CSV
    df.to_csv(args.output, float_format="%.6f")
    print(f"\nSaved table to {args.output}")

    ### how to extract x and y for plotting
    x = df.loc["distance"].to_numpy(dtype=float)              # (m,)

    ### selected rows
    #selected_rows = ["Kinetic","Total"]
    #y_total = df.loc["Total"].to_numpy(dtype=float)   # (m,)
    if args.selected_rows[0].isdigit():
        selected_rows = []
        for index in args.selected_rows:
            selected_rows.append(ENERGY_KEYS[int(index)])
    else:
        selected_rows = args.selected_rows
    y_all = df.loc[selected_rows].to_numpy()    # (n, m)

    plot_dict = {}
    plot_dict['legends'] = selected_rows
    mplot_nvector(x, y_all, plot_dict=plot_dict)

    return 0

if __name__ == "__main__":
    main()
