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
    - run /qcfs/noische/scripts/test/lammps_gethbonds.py for pkl_file
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
    hbonds = pickle.load(f)

print("loading pickle %s .." % layer_file)
with open(layer_file) as f:
    layers = pickle.load(f)

timesteps = sorted(hbonds.keys())

result = dict()

def find_layer(z):
    if z_bin:
        return int(np.digitize(z, z_bin))
    else:
        return 0

# bin
for t in tqdm.tqdm(timesteps, ncols=120, desc="Hbonds"):
    # count natoms in each layer
    layer_id = list(set(layers[t].values()))
    n_atoms = [0 for i in layer_id]
    for i in layers[t]:
        n_atoms[layers[t][i]] += 1

    hb = dict() # stores hbonds at time t
    avg_hb = [0 for i in layer_id] # average hbonds inside the layer

    for i in hbonds[t]:
        for j in i[:1]:
            layer = (find_layer(i[2]), find_layer(i[3]))
            if layer in hb:
                hb[layer] += 1
            else:
                hb[layer] = 1

    # average number of hbonds "inside" the layer (intralayer)
    for i in layer_id:
        avg_hb[i] = hb[(i, i)]/float(n_atoms[i])

    result[t] = avg_hb

# write output
output = "#t " + " ".join(str(i) for i in layer_id) + "\n"
for t in timesteps:
    output += "%d " % t + " ".join(str(i) for i in result[t]) + "\n"

with open("hbonds.layer.dat", 'w') as f:
    f.write(output)

