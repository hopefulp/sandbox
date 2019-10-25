#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys
import pickle
import numpy as np
import tqdm
import bgf
import nutils as nu

# command parameters
"""
hbonds[t] = [d_atom.aNo, a_atom.aNo, [donor.x, y, z], [acceptor.x, y, z], dist, theta]
"""
usage = """%s bgf_file pkl_file (selection) (layer_file)
    prerequisites:
    - run /qcfs/noische/scripts/test/lammps_gethbonds.py for pkl_file
    - run /qcfs/noische/scripts/test/lammps_count_layer.py for layer_file
""" % sys.argv[0]

if len(sys.argv) < 2:
    print(usage)
    sys.exit(0)

bgf_file = sys.argv[1]
pkl_file = sys.argv[2]

if len(sys.argv) == 3:
    selection = ""
    layer_file = ""
elif len(sys.argv) == 4:
    selection = ""
    layer_file = sys.argv[3]
else:
    selection = sys.argv[3]
    layer_file = sys.argv[4]


# init
mybgf = bgf.BgfFile(bgf_file)
dist_delta = 0.05
angle_delta = 0.2
Ehb_delta = 0.2
bin_dist = np.arange(2, 4.0 + dist_delta, dist_delta)
bin_angle = np.arange(0, 30.0 + angle_delta, angle_delta)
bin_Ehb = np.arange(-15, 15.0 + Ehb_delta, Ehb_delta)

#alpha_delta = 1.0   # v4
#bin_alpha = np.arange(-90.0, 90.0 + alpha_delta, alpha_delta)  # v4


# pickle load
print("loading pickle %s .." % pkl_file)
with open(pkl_file, 'rb') as f:
    hbonds = pickle.load(f)

timesteps = sorted(hbonds.keys())

layer = dict(); layer_id = [0]
if layer_file:
    print("reading layer file %s .." % layer_file)
    with open(layer_file, 'rb') as f:
        layer = pickle.load(f)
    
    for t in timesteps:
        try:
            layer_id = list(set(layer[t].values()))
        except KeyError:
            nu.die('Timestep t not found in the layer file.')

### Variables
# per t
t_n_hbonds = dict();    # total number of HBonds
t_n_atoms = dict();     # total number of atoms in the selection
t_hb_atoms = dict();    # total number of atoms which have HBonds
t_avg_Ehb = dict();     # average E_hbond value ##

t_n_atom_layer = dict();    # number of atoms in a layer
t_n_hb_layer = dict();      # number of HBonds in a layer 
t_Ehb_layer = dict();   # average E_hbond value per layer ##

t_distr_dist = dict();  # distributions of HBond distance
t_distr_angle = dict(); # distributions of HBond angle
t_distr_Ehb = dict();    # distributions of E_hbond ##
#t_distr_alpha = dict(); # distributions of angle between H-O-H and O..O # v4

# average over t
avg_n_atom_layer = dict();  # average number of atoms in each layer
avg_n_hb_atom_layer = dict();   # average number of HBonds?????

avg_distr_dist = [];    # average distribution of HBond distances
avg_distr_angle = [];   # average distribution of HBond angles
avg_distr_Ehb = [];     # average distribution of E_hbond ##
#avg_distr_alpha = [];   # average distribution of H-O-H and O..O angles # v4

### looping over t
for t in tqdm.tqdm(timesteps, ncols=120, desc="Hbonds"):
    ### Init
    anos = []
    if layer_file:
        anos = layer[t].keys()
    else:
        temp = dict()
        for atom in mybgf.a:
            if selection:
                if "O" in atom.ffType and eval(selection):
                    anos.append(atom.aNo)
            else:
                if "O" in atom.ffType:
                    anos.append(atom.aNo)
        for ano in anos:
            temp[ano] = 0
        layer[t] = temp
    n_atoms = len(anos)

    # per atom
    n_hb = dict()
    for ano in anos:
        n_hb[ano] = 0
    n_hbonds = 0
    hb_atoms = 0

    # per system
    dist = []
    angle = []
    Ehb = []
    #alpha = []  # H-O-H and O..O angle  # v4

    # per layer
    n_atom_layer = dict()  
    n_hb_layer = dict()
    Ehb_layer = dict()
    for i in layer_id:
        n_atom_layer[i] = 0
        n_hb_layer[i] = 0
        Ehb_layer[i] = 0.0

    # assign layer id to atom
    for i in hbonds[t]:
        n_hb[i[0]] += 1
        n_hb[i[1]] += 1
        dist.append(i[4])
        angle.append(i[5])
        Ehb.append(i[6][2])
        #alpha.append(i[8])

    distr_dist, _ = np.histogram(dist, bins=bin_dist)
    distr_angle, __ = np.histogram(angle, bins=bin_angle)
    distr_Ehb, ___ = np.histogram(Ehb, bins=bin_Ehb)
    #distr_alpha, ____ = np.histogram(alpha, bins=bin_alpha)
    norm_distr_dist = distr_dist/float(distr_dist.sum())
    norm_distr_angle = distr_angle/float(distr_angle.sum())
    norm_distr_Ehb = distr_Ehb/float(distr_Ehb.sum())
    #norm_distr_alpha = distr_alpha/float(distr_alpha.sum())

    # count hbonds and natoms per layer
    for ano in anos:
        n_hbonds += n_hb[ano]
        if n_hb[ano] != 0: hb_atoms += 1
        n_atom_layer[layer[t][ano]] += 1
        n_hb_layer[layer[t][ano]] += n_hb[ano]

    # store
    t_n_atoms[t] = n_atoms
    t_n_hbonds[t] = n_hbonds
    t_hb_atoms[t] = hb_atoms
    t_distr_dist[t] = distr_dist
    t_distr_angle[t] = distr_angle
    t_n_atom_layer[t] = n_atom_layer
    t_n_hb_layer[t] = n_hb_layer
    t_avg_Ehb[t] = np.mean(np.array(Ehb))

    avg_distr_dist.append(norm_distr_dist)
    avg_distr_angle.append(norm_distr_angle)
    avg_distr_Ehb.append(norm_distr_Ehb)
    #avg_distr_alpha.append(norm_distr_alpha)

# Average over t
avg_distr_dist = np.mean(avg_distr_dist, axis=0)
avg_distr_angle = np.mean(avg_distr_angle, axis=0)
avg_distr_Ehb = np.mean(avg_distr_Ehb, axis=0)
#avg_distr_alpha = np.mean(avg_distr_alpha, axis=0)

# Write output
with open(pkl_file + ".hbonds.analysis.dat", 'w') as f:
    output = "#t\t"
    output += "total hb\t"
    output += "total atoms\t"
    output += "hbd atoms\t"
    output += "nhb/atoms\t"
    output += "".join("nhb/layer %d\t" % int(i+1) for i in layer_id)
    output += "E_hbond" + "\n"

    for t in timesteps:
        output += "%d\t" % t
        output += "%d\t" % t_n_hbonds[t]
        output += "%d\t" % t_n_atoms[t]
        output += "%d\t" % t_hb_atoms[t]
        output += "%8.3f\t" % (float(t_n_hbonds[t])/float(t_n_atoms[t]))
        output += "".join("%8.3f\t" % (float(t_n_hb_layer[t][i])/float(t_n_atom_layer[t][i])) for i in layer_id)
        output += "%8.3f" % t_avg_Ehb[t] + "\n"

    f.write(output)

with open(pkl_file + ".hbonds.analysis.dist.dat", 'w') as f:
    output = "#d(O-O)\tfrequency\t#delta: %s\n" % str(dist_delta)
    #for i, j in zip(bin_dist[1:], avg_distr_dist):
    for i, j in zip(_[1:], avg_distr_dist):
        output += "%5.3f" % (i-dist_delta/2.0) + "%8.5f" % j + "\n"
    f.write(output)

with open(pkl_file + ".hbonds.analysis.theta.dat", 'w') as f:
    output = "#theta(H-O-O)\tfrequency\t#delta: %s\n" % str(angle_delta)
    #for i in zip(bin_angle[1:], avg_distr_angle):
    for i, j in zip(__[1:], avg_distr_angle):
        output += "%5.3f" % (i-angle_delta/2.0) + "%8.5f" % j + "\n"
    f.write(output)

with open(pkl_file + ".hbonds.analysis.E_hbond.dat", 'w') as f:
    output = "#E_hbond\tfrequency\t#delta: %s\n" % str(Ehb_delta)
    #for i in zip(bin_Ehb[1:], avg_distr_Ehb):
    for i, j in zip(___[1:], avg_distr_Ehb):
        output += "%5.3f" % (i-Ehb_delta/2.0) + "%8.5f" % j + "\n"
    f.write(output)

'''
with open(pkl_file + ".hbonds.analysis.alpha.dat", 'w') as f:
    output = "#alpha(H-O-H - O..O)\tfrequency\t#alpha: %s\n" % str(alpha_delta)
    #for i in zip(bin_alpha[1:], avg_distr_alpha):
    for i, j in zip(____[1:], avg_distr_alpha):
        output += "%5.3f" % (i-alpha_delta/2.0) + "%8.5f" % j + "\n"
    f.write(output)
'''


