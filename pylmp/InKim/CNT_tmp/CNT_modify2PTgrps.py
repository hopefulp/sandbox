#!/home/noische/python

import sys
import os
import time

from numpy import pi
import bgf
import CNT
from LAMMPS_trj2bgf import *

usage = """
CNT_modify2PTgrps.py bgf_file grps_file 

    This script modifies GroupVolume part in grps_file.
    This script should be run within Tod's LAMMPS-2PT pbs script.
"""
version = "150803"

if len(sys.argv) < 5:
    print usage;
    sys.exit(0);

bgf_filename = sys.argv[1]
grps_filename = sys.argv[2]

# read LAMMPS trajectory from scratch
nt = CNT.Nanotube(bgf_filename)
nt.calc_height_radius() # calculate effective volume within the nanotube
if nt.type == "BNT":
    vdw_radii = (1.92 + 1.55) / 2    # B: http://periodic.lanl.gov/5.shtml, N: http://periodic.lanl.gov/7.shtml
    vdw_radii = (1.790706425 + 1.6312801) / 2    # B and N: lammps data file
else:
    vdw_radii = 3.38383824 / 2    # C_2G in QMFF# CNT
    vdw_radii = 2.45    # C_2G from trajectory
height = nt.height + 2 * vdw_radii
total_volume = nt.pbc[0] * nt.pbc[1] * nt.pbc[2]

inner_volume = pi * (nt.radius - vdw_radii)**2 * height    # pi * r**2 * h
nt_volume = pi * ( (nt.radius + vdw_radii)**2 - (nt.radius - vdw_radii)**2 ) * height    # pi * (ro**2 - ri**2) * h
outer_volume = total_volume - inner_volume - nt_volume

# write GroupVolume to grps file
g = grpfile(grps_filename)
g.grp[1]['volume'] = nt_volume
g.grp[2]['volume'] = inner_volume
g.grp[3]['volume'] = outer_volume
g.write(grps_filename, zip=True)

grps_dir = os.path.dirname(grps_filename)
logf = open(grps_dir + "/" + "CNT_modify2PTvol.log", 'a')
logf.write(str(time.asctime(time.gmtime())) + '\t' + str(timestep) + '\t' + str(nt.radius) + '\t' + str(vdw_radii) + '\t' + str(yes_exterior_water) + '\t' + str(total_volume) + '\t' + str(nt_volume) + '\t' + str(inner_volume) + '\t' + str(outer_volume) + '\n')
logf.close()

### end of code
