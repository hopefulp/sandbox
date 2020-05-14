#!/usr/bin/env python3

import sys

import matplotlib.pyplot as plt

import lammps_log as llog

usage = """
usage: lammps_log_plot.py log_file thermo_keyword"""

if len(sys.argv) < 2:
    print(usage)
    sys.exit(1)

log_file = sys.argv[1]
thermo_keyword = sys.argv[2:]

mylog = llog.Log(log_file)
mylog.plot(thermo_keyword)

