#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys
import cPickle as pickle
import numpy as np
import tqdm
import nutils as nu

"""
aop[t] = [angle, i.z, atom.z, j.z]
"""
usage = """%s pkl_file 'z_bin_ranges'
    suggestions: (for Graphene/MoS2 analysis)
    
    - graphene single layer: z_bin_ranges = ""
    - graphene double layer: z_bin_ranges = "20"
    - graphene triple layer: z_bin_ranges = "18.7 21.2"
""" % sys.argv[0]

pkl_file = sys.argv[1]
z_range = sys.argv[2]

# init
angle_bin = np.arange(0, 181, 1.0)
z_bin = [float(i) for i in z_range.split()]

# load
print("loading pickle %s .." % pkl_file)
with open(pkl_file) as f:
    aop = pickle.load(f)

timesteps = sorted(aop.keys())
avg_total_angles = [];
avg_intra_angles = [];
avg_inter_angles = [];
avg_center_angles = [];

# bin
for t in tqdm.tqdm(timesteps, ncols=120, desc="Binning"):
    total_angles = []
    intra_angles = []
    inter_angles = []
    center_angles = []
    for i in aop[t]:
        total_angles.append(i[0])
        if len(z_bin) != 0:
            if np.digitize(i[1], z_bin) == np.digitize(i[2], z_bin) and np.digitize(i[2], z_bin) == np.digitize(i[3], z_bin):
                # intra layer
                intra_angles.append(i[0])
            else:
                # inter layer
                inter_angles.append(i[0])
                if len(z_bin) > 1 and np.digitize(i[2], z_bin) == 1:
                    # center layer
                    center_angles.append(i[0])

    hist_total_angle, a = np.histogram(total_angles, bins=angle_bin)
    hist_intra_angle, a = np.histogram(intra_angles, bins=angle_bin)
    hist_inter_angle, a = np.histogram(inter_angles, bins=angle_bin)
    hist_center_angle, a = np.histogram(center_angles, bins=angle_bin)
    avg_total_angles.append(hist_total_angle/float(len(total_angles)))
    avg_intra_angles.append(hist_intra_angle/float(len(total_angles)))
    avg_inter_angles.append(hist_inter_angle/float(len(total_angles)))
    avg_center_angles.append(hist_center_angle/float(len(total_angles)))

# average
#avg_angles = np.mean(np.array(avg_angles), axis=0)
avg_total_angles = np.mean(np.array(avg_total_angles), axis=0)
avg_intra_angles = np.mean(np.array(avg_intra_angles), axis=0)
avg_inter_angles = np.mean(np.array(avg_inter_angles), axis=0)
avg_center_angles = np.mean(np.array(avg_center_angles), axis=0)

# print
print("#degree total intra inter")
for i in zip(a, avg_total_angles, avg_intra_angles, avg_inter_angles, avg_center_angles):
    print(" ".join(str(x) for x in i))
