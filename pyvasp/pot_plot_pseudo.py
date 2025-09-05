#!/home/joonho/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt

outf = open('fft.dat', 'w')

# Load V(G) from POTCAR (real, only positive Gs)
v_g_pos = np.loadtxt('pseudo_k.dat')  # Or from your uploaded array
v_g = v_g_pos[:, 1]  # assuming column 1 has local potential

# Reconstruct Hermitian-symmetric spectrum (positive + negative Gs)
v_g_sym = np.concatenate((v_g[-2:0:-1], v_g[:] ))  # Exclude the first and last to avoid double counting

# Apply inverse FFT
#v_r = np.fft.ifftshift(v_g_sym) # Ensure real result
#v_r = np.fft.ifft(np.fft.ifftshift(v_g_sym)).real # Ensure real result
v_r = np.fft.ifft(np.fft.ifftshift(v_g_sym)) # Ensure real result
v_r = np.real(v_r)
#v_r = np.fft.fftshift(v_r).real
# Define real-space grid
N = len(v_r)
a = 1.0  # total length of real-space domain (arbitrary unit unless you scale with Gmax)
r = np.linspace(0, a, N, endpoint=False)

for x, y in zip(r, v_r):
    outf.write(f"{x:>15.10f} {y:>15.10f} \n")

# Plot
plt.plot(r, v_r)
plt.title("Reconstructed Real-space Local Potential from POTCAR V(G)")
plt.xlabel("r (arbitrary units)")
plt.ylabel("V(r)")
plt.grid(True)
plt.tight_layout()
plt.show()
