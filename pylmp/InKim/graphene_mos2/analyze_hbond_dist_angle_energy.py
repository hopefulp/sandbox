#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys
import os
import cPickle as pickle
import numpy as np
import tqdm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.mlab import griddata

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

test = False 

# init
dr = 0.02; db = 1
if test: 
    dr = 0.1; db = 5.0

bin_dist = np.arange(2.5, 3.5, dr)
bin_angle = np.arange(0 + db, 30, db)
dist = []
angle = []
E_hbond = []

if not os.path.exists("hbond.dist-angle-Ehbond.pickle"):
    # pickle load
    print("loading pickle %s .." % pkl_file)
    with open(pkl_file) as f:
        hbonds = pickle.load(f)

    timesteps = sorted(hbonds.keys())

    ### looping over t
    for t in tqdm.tqdm(timesteps, ncols=120, desc="Hbonds"):
        # assign layer id to atom
        for i in hbonds[t]:
            dist.append(i[4])
            angle.append(i[5])
            E_hbond.append(i[6][2])

    # store
    with open("hbond.dist-angle-Ehbond.pickle", 'wb') as f:
        print("writing pickle hbond.dist.pickle")
        pickle.dump([dist, angle, E_hbond], f)
else:
    with open("hbond.dist-angle-Ehbond.pickle", 'rb') as f:
        print("loading pickle hbond.dist.pickle")
        dist, angle, E_hbond = pickle.load(f)

print('Gridding..')
Z = griddata(dist, angle, E_hbond, bin_dist, bin_angle, interp="linear")
X, Y = np.meshgrid(bin_dist, bin_angle)


print('Plotting..')
bounds = np.array([-8.0, -7.0, -6.0, -5.0, -4.0, -3.0])
CS = plt.contourf(X, Y, Z, bounds)
cbar = plt.colorbar(CS)
cbar.ax.set_ylabel('E_hbond')
plt.show()


'''
fig = plt.figure()
H, xedges, yedges = np.histogram2d(dist, angle, bins=(bin_dist, bin_angle), normed=True)
H = H.transpose()

r, b = np.meshgrid(xedges[1:], yedges[1:])

rho = 55.5 * 6.022 * 0.0001 # water number density
norm = rho * 2.0 * np.pi * np.sin(np.radians(b)) * np.radians(db) * r**2 * dr
g_RB = 10.0 * H / norm / 658    # why??
pmf_RB = -np.log(g_RB)
pmf_RB_mask = np.ma.masked_where(np.isinf(pmf_RB), pmf_RB)
if test: 
    print("%s\n\n" % str(H))
    print("%s\n\n" % str(r))
    print("%s\n\n" % str(b))
    print("%s\n\n" % str(norm))
    print("%s\n\n" % str(g_RB))
    print("%s\n\n" % str(pmf_RB_mask))

min, max = -np.abs(pmf_RB_mask).max(), np.abs(pmf_RB_mask).max()
min, max = -3.0, 3.0
bounds = np.array([-3.0, -2.0, -1.0, 0, 1.0, 1.23, 2.0, 3.0])
contourcolor = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
plt.pcolormesh(xedges, yedges, pmf_RB_mask, cmap='RdBu_r', vmin=min, vmax=max, norm=contourcolor)
#plt.pcolormesh(xedges, yedges, pmf_RB, cmap='RdBu_r', vmin=min, vmax=max, norm=contourcolor)

plt.colorbar()
plt.show()
'''
