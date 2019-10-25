#!/qcfs/noische/program/psi/python/3.5.2/bin/python3.5
import sys
import numpy as np
import pandas as pd
import pickle
import tqdm
import nutils as nu

usage = """This script is used to analyze H-bond calculation results from gethbond_parallel.py (v5)
Usage: %s pickle_file
""" % sys.argv[0]
# pickle file
pkl_file = sys.argv[1]

# histogram
dist_delta = 0.01
angle_delta = 0.05
Ehb_delta = 0.05

bin_dist = np.arange(2, 4.0 + dist_delta, dist_delta)
bin_angle = np.arange(0, 30.0 + angle_delta, angle_delta)
bin_Ehb = np.arange(-15, 15.0 + Ehb_delta, Ehb_delta)

# data load
with open(pkl_file, 'rb') as f:
    data = pickle.load(f)

ts = sorted(data.keys())

avg_nhb = []
avg_distr_dist = []
avg_distr_angle = []
avg_distr_Ehb = []
avg_dist = []
avg_angle = []
avg_Ehb = []
for t in tqdm.tqdm(ts, ncols=120):
    nhb = dict()
    natoms = len(data[t][0])
    df = pd.DataFrame(data[t][1])

    # n_hbonds (0 & 1)
    #df = df.drop_duplicates(subset=4, keep='first')
    #nhb = [len(df[df[0] == i]) + len(df[df[1] == i]) for i in df[0]]
    #nhb = [len(df[df[0] == i]) for i in df[0]]

    hb = dict()
    for i in zip(df[0], df[1]):
        for j in i:
            if j in hb.keys():
                hb[j] += 1
            else:
                hb[j] = 1

    nhb = [hb[i] for i in hb.keys() if hb[i]]

    nhb = np.array(nhb)
    #avg_nhb.append([natoms, nhb.sum(), nhb.sum()/natoms])
    avg_nhb.append([natoms, nhb.sum(), nhb.mean()])

    # dist (2)
    distr_dist, _ = np.histogram(df[2], bins=bin_dist)
    norm_distr_dist = distr_dist/float(distr_dist.sum())
    avg_distr_dist.append(norm_distr_dist)
    avg_dist.append(df[2].mean())

    # angle (3)
    distr_angle, __ = np.histogram(df[3], bins=bin_angle)
    norm_distr_angle = distr_angle/float(distr_angle.sum())
    avg_distr_angle.append(norm_distr_angle)
    avg_angle.append(df[3].mean())

    # E_hbond dist (4)
    distr_Ehb, ___ = np.histogram(df[4], bins=bin_Ehb)
    norm_distr_Ehb = distr_Ehb/float(distr_Ehb.sum())
    avg_distr_Ehb.append(norm_distr_Ehb)
    avg_Ehb.append(df[4].mean())

#avg_nhb = np.mean(avg_nhb, axis=0)
avg_distr_dist = np.mean(avg_distr_dist, axis=0)
avg_distr_angle = np.mean(avg_distr_angle, axis=0)
avg_distr_Ehb = np.mean(avg_distr_Ehb, axis=0)

avg_dist = np.mean(avg_dist)
avg_angle = np.mean(avg_angle)
avg_Ehb = np.mean(avg_Ehb)

with open(pkl_file + ".hbonds.analysis.average.dat", 'w') as f:
    output = "natoms\tsum(hb)\thb_per_atom\n"
    for i in avg_nhb:
        output += "%d %d %f\n" % (i[0], i[1], i[2])
    f.write(output)

with open(pkl_file + ".hbonds.analysis.dist.dat", 'w') as f:
    output = "d\tfrequency\n"
    for i, j in zip(_[1:], avg_distr_dist):
        output += "%5.3f" % (i-dist_delta/2.0) + "%8.5f" % j + "\n"
    f.write(output)

with open(pkl_file + ".hbonds.analysis.theta.dat", 'w') as f:
    output = "theta\tfrequency\n"
    for i, j in zip(__[1:], avg_distr_angle):
        output += "%5.3f" % (i-angle_delta/2.0) + "%8.5f" % j + "\n"
    f.write(output)

with open(pkl_file + ".hbonds.analysis.E_hbond.dat", 'w') as f:
    output = "E_hbond\tfrequency\n"
    for i, j in zip(___[1:], avg_distr_Ehb):
        output += "%5.3f" % (i-Ehb_delta/2.0) + "%8.5f" % j + "\n"
    f.write(output)

avg_nhb = np.mean(avg_nhb, axis=0)


print("Average atoms: {0: 8.3f} Total H-bonds: {1:8.3f} Average H-bond numbers: {2:8.3f} / molecule".format(avg_nhb[0], avg_nhb[1], avg_nhb[2]))
print('Average H-bond distance: {0:8.3f}'.format(avg_dist))
print('Average H-bond angle: {0:8.3f}'.format(avg_angle))
print('Average H-bond energy: {0:8.3f}'.format(avg_Ehb))
