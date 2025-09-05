#!/home/joonho/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def load_potential(filename):
    r = []
    V = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    r.append(float(parts[0]))
                    V.append(float(parts[1]))
                except ValueError:
                    continue  # skip lines that can't be parsed
    return np.array(r), np.array(V)


# Load potentials
r1, V_ae = load_potential("VTOTAL.ae")
r2, V_occ = load_potential("VTOTAL_OCC")


# Load potentials
r1, V_ae = load_potential("VTOTAL.ae")
r2, V_occ = load_potential("VTOTAL_OCC")

# Interpolate VTOTAL_OCC onto r1 grid
interp_func = interp1d(r2, V_occ, bounds_error=False, fill_value="extrapolate")
V_occ_interp = interp_func(r1)

# Now you can compute the difference safely
V_diff = V_ae - V_occ_interp

# Plotting
plt.figure(figsize=(10, 6))
#plt.plot(r1, V_ae, label='VTOTAL.ae (Neutral)', color='blue')
#plt.plot(r2, V_occ, label='VTOTAL_OCC (Half-Ionized)', color='green')
plt.plot(r1, V_diff, label='Difference (ae - occ)', color='red', linestyle='--')

plt.xlabel('Radius r (Bohr)')
plt.ylabel('Potential (Hartree)')
plt.title('Comparison of Atomic Potentials')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("potential_comparison.png", dpi=300)
plt.show()

