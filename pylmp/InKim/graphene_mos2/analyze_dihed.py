#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys
import cPickle as pickle
import numpy as np
import tqdm
import nutils as nu

"""
hbonds[t] = [d_atom.aNo, a_atom.aNo, d_atom.z, a_atom.z]
"""
usage = """%s pkl_file layer_file 'z_bin_ranges'
    single layer: z_bin_ranges = ""
    double layer: z_bin_ranges = "20" for graphene
    triple layer: z_bin_ranges = "18.7 21.2" for graphene

    prerequisites:
    - run /qcfs/noische/scripts/test/lammps_get_aop.py for aop pkl_file
    - run /qcfs/noische/scripts/test/lammps_count_layer.py for layer_file
""" % sys.argv[0]

if len(sys.argv) < 3:
    print(usage)
    sys.exit(0)

pkl_file = sys.argv[1]
layer_file = sys.argv[2]
z_range = sys.argv[3]

# init
angle_bin = np.arange(0, 181, 1.0)
z_bin = [float(i) for i in z_range.split()]

# load
print("loading pickle %s .." % pkl_file)
with open(pkl_file) as f:
    diheds = pickle.load(f)

print("loading pickle %s .." % layer_file)
with open(layer_file) as f:
    layers = pickle.load(f)

timesteps = sorted(diheds.keys())

result = dict()

def find_layer(z):
    if z_bin:
        return int(np.digitize(z, z_bin))
    else:
        return 0

avg_hist_dihed = dict()

# bin
for t in tqdm.tqdm(timesteps, ncols=120, desc="Hbonds"):
    # count natoms in each layer
    layer_id = list(set(layers[t].values()))
    n_atoms = [0 for i in layer_id]
    if not n_atoms:
        n_atoms = [0]
    if not layer_id:
        layer_id = [0]
    for i in layers[t]:
        n_atoms[layers[t][i]] += 1

    d = dict() # stores diheds at time t
    for i in diheds[t]:
        layer = (find_layer(i[1]), find_layer(i[2]))
        if layer in d:
            d[layer].append(i[0])
        else:
            d[layer] = [i[0]]

    # d[(0,0)] = dihed angles in layer 0 // d[(0, 1)] = dihed angles lying across 0 and 1 // ...
    # average number of diheds "inside" the layer (intralayer)
    for i in layer_id:
        hist_dihed, _ = np.histogram(d[(i, i)], bins=angle_bin, density=True)
        if (i, i) in avg_hist_dihed:
            avg_hist_dihed[(i, i)].append(hist_dihed)
        else:
            avg_hist_dihed[(i, i)] = [hist_dihed]

for i in avg_hist_dihed:
    avg_hist_dihed[i] = np.mean(np.array(avg_hist_dihed[i]), axis=0)

# write output
output = ""
for i in avg_hist_dihed:
    output += "#angle layer%s\n" % str(i)
    for j in zip(_, avg_hist_dihed[i]):
        output += " ".join(str(k) for k in j) + "\n"
    output += "\n"

with open("diheds.layer.dat", 'w') as f:
    f.write(output)

