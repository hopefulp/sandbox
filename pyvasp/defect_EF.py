#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# Fermi level range (0 → band gap)
EF = np.linspace(0, 2.0, 200)

# Charge states: formation energy at EF=0 and slopes
charge_states = {
    0: {'E0': 0.5, 'q': 0},
    1: {'E0': 1.0, 'q': 1},
    2: {'E0': 1.8, 'q': 2}
}

plt.figure(figsize=(6,4))

for q, data in charge_states.items():
    E_f = data['E0'] + data['q']*EF
    plt.plot(EF, E_f, label=f'q={q}')

# Mark transitions (intersections)
# ε(0/+1)
epsilon_01 = (charge_states[1]['E0'] - charge_states[0]['E0']) / (charge_states[0]['q'] - charge_states[1]['q'])
# ε(+1/+2)
epsilon_12 = (charge_states[2]['E0'] - charge_states[1]['E0']) / (charge_states[1]['q'] - charge_states[2]['q'])

plt.axvline(epsilon_01, color='gray', linestyle='--', label='ε(0/+1)')
plt.axvline(epsilon_12, color='gray', linestyle=':', label='ε(+1/+2)')

plt.xlabel('Fermi level E_F (eV, relative to VBM)')
plt.ylabel('Defect formation energy E_f (eV)')
plt.title('Defect Formation Energy vs Fermi Level')
plt.legend()
plt.grid(True)
plt.show()
