#!/home/noische/python
"""
LAMMPS_density.py -b bgf_file -t trj_file -f ff_file
read bgf -> get mass
read trj -> get volume
calculates system density
"""

import sys
import functools
import getopt
import pandas as pd
import bgf
import bgftools as bt
import nutils as nu
#import lammpstrj as lt
from lammps.trj import *

version = "161130"

def get_volume(trj_file):
    #mytrj = lt.lammpstrj(trj_file)
    mytrj = Trj(trj_file)
    mytrj.load()
    ts = sorted(mytrj.timesteps)
    return {t: functools.reduce(lambda x, y: x * y, mytrj.pbc[t]) for t in ts}


def get_mass(bgf_file, ff_file=''):
    mybgf = bgf.BgfFile(bgf_file)

    return bt.getMass(mybgf, ff_file=ff_file)


def get_density(bgf_file, trj_file, ff_file="", draw=False, out_file="", n=0, silent=False):
    mass = get_mass(bgf_file, ff_file=ff_file)  # float
    volume = get_volume(trj_file)   # dict
    density = {d: mass/volume[d]/6.022*10 for d in volume.keys()}

    df = pd.DataFrame(density, index=['density'])
    df = df.T

    if n < len(df):
        df = df[-n:]
    else:
        nu.warn("Not enough sample shots to calculate density with requested number of last snapshots.")
    
    density = float(df.mean())

    if not silent:
        print(df)
        print("Average density: {0:<11.5f}".format(density))

    if draw:
        import matplotlib.pyplot as plt
        plt.figure(); df.plot(); plt.legend(loc='best'); fr = plt.gca(); fr.axes.set_ylim(bottom=0);
        plt.show()

    if out_file:
        df.to_csv(out_file, sep=' ', index_label='t')
        nu.warn("the character # should be added to the first of the line 1 if you want to plot the file on gnuplot.")

    return density


if __name__ == "__main__":

    bgf_file = ""; trj_file = ""; ff_file = ""; flag_draw = False; out_file = ""; nsamples = 0

    usage = """
%s -b bgf_file -t trj_file [OPTIONS]
    Calculates density of the system in a LAMMPS trajectory file.

    Additional options:
        -f ff_file      force field files to calculate accurate atomic masses
        -o out_file     a filename to save the result in text format
        -d              draws a plot
        -n              specify the number of last snapshots to average 

    Please report any bugs to in.kim@kaist.ac.kr.
    Last updated: %s
        """ % (sys.argv[0], version)

    if len(sys.argv) < 2:
        print(usage)
        sys.exit(0)

    options, args = getopt.getopt(sys.argv[1:], 'hb:t:f:do:n:', ['help', 'bgf=', 'trj=', 'ff=', 'draw=', 'out=', 'n='])

    print(sys.argv)
    print("Requested options: " + str(options))

    for option, value in options:
        if option in ('-h', '--help'):
            print(usage)
            sys.exit(0)
        elif option in ('-b', '--bgf'):
            bgf_file = value
        elif option in ('-t', '--trj'):
            trj_file = value
        elif option in ('-f', '--ff'):
            ff_file = str(value).strip()
        elif option in ('-d', '--draw'):
            flag_draw = True
        elif option in ('-o', '--out'):
            out_file = str(value)
        elif option in ('-n', '--n'):
            nsamples = int(value)
        elif option == NULL:
            print(usage)
            sys.exit(0)

    get_density(bgf_file, trj_file, ff_file=ff_file, draw=flag_draw, out_file=out_file, n=nsamples)
