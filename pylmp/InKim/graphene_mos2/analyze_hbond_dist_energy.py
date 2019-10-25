#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys
import os
import cPickle as pickle
import numpy as np
import tqdm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import NullFormatter

# command parameters
"""
hbonds[t] = [d_atom.aNo, a_atom.aNo, [donor.x, y, z], [acceptor.x, y, z], dist, theta]
"""
usage = """%s bgf_file pkl_file (layer_file)
    prerequisites:
    - run /qcfs/noische/scripts/test/lammps_gethbonds.py for pkl_file
    - run /qcfs/noische/scripts/test/lammps_count_layer.py for layer_file
""" % sys.argv[0]

if len(sys.argv) < 2:
    print(usage)
    sys.exit(0)

pkl_file = sys.argv[1]


# init
dist = []
E_hbond = []

# pickle load
if not os.path.exists("hbond.dist-E_hbond.pickle"):
    print("loading HBond data pickle %s .." % pkl_file)
    with open(pkl_file) as f:
        hbonds = pickle.load(f)

    timesteps = sorted(hbonds.keys())

    ### Variables

    ### looping over t
    for t in tqdm.tqdm(timesteps, ncols=120, desc="Hbonds"):
        # assign layer id to atom
        for i in hbonds[t]:
            dist.append(i[4])
            E_hbond.append(i[6][2])

    # store
    with open("hbond.dist-E_hbond.pickle", 'wb') as f:
        print("writing pickle hbond.dist-E_hbond.pickle")
        pickle.dump([dist, E_hbond], f)
else:
    with open("hbond.dist-E_hbond.pickle", 'rb') as f:
        print("Reading pickle hbond.dist-E_hbond.pickle")
        dist, E_hbond = pickle.load(f)

# matplotlib decoration
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left + width + 0.02

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]

# rectangular figure
plt.figure(1, figsize=(8, 8))
axScatter = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx)
axHisty = plt.axes(rect_histy)
axScatter.scatter(dist, E_hbond, marker='.')
axScatter.set_xlabel('O-O dist')
axScatter.set_ylabel('E_hbond')

# histograms
binwidth_dist = 0.02
binwidth_E = 0.5
bin_dist = np.arange(2.4, 3.5, binwidth_dist); dist_min, dist_max = np.min(bin_dist), np.max(bin_dist)
bin_E = np.arange(-15, 10, binwidth_E); E_min, E_max = np.min(bin_E), np.max(bin_E)
axScatter.set_xlim((dist_min, dist_max))
axScatter.set_ylim((E_min, E_max))

axHistx.hist(dist, bins=bin_dist, normed=1)
axHisty.hist(E_hbond, bins=bin_E, orientation='horizontal', normed=1)
axHistx.set_xlim(axScatter.get_xlim())
axHisty.set_ylim(axScatter.get_ylim())

#plt.plot(dist, E_hbond, '.')
plt.show()
