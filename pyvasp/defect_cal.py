import numpy as np
import matplotlib.pyplot as plt

# --- 입력(DFT 계산값: 예시값) ---
E_bulk = -102.5           # total energy of pristine supercell (eV)
E_def_q = {0: -102.0, 1: -100.5, 2: -99.0}  # total energy of defect supercell for q=0,1,2
E_O2 = -9.86             # DFT energy of O2 molecule (eV)
E_AxOy = -110.0          # DFT energy of A_xO_y compound (e.g., AO)
E_A_ref = -50.0          # DFT energy per atom of pure A (reference)

x=1; y=1   # stoichiometry AO (example)

# --- compute mu_O range ---
muO_ref = 0.5 * E_O2   # reference
# O-rich: delta_muO = 0 -> muO = muO_ref
muO_Orich = muO_ref

# O-poor lower bound from mu_A = mu_A_ref:
delta_muO_min = (E_AxOy - x*E_A_ref)/y - muO_ref
muO_opoor = muO_ref + delta_muO_min

print("muO_ref, O-rich, O-poor:", muO_ref, muO_Orich, muO_opoor)

# choose one condition
muO = muO_Orich   # or muO_opoor

# compute mu_A from compound
muA = (E_AxOy - y*muO)/x

# --- compute formation energies vs E_F ---
Ev = 0.0
Eg = 2.0  # band gap (if needed)
EFs = np.linspace(0, Eg, 201)
Ef_vs_EF = {q: [] for q in E_def_q.keys()}

# assume defect removes one O: n_O = -1
n_O = -1

for EF in EFs:
    for q, Edef in E_def_q.items():
        Ef = Edef - E_bulk - n_O*muO + q*(EF + Ev)  # + E_corr omitted
        Ef_vs_EF[q].append(Ef)

# plot
for q in Ef_vs_EF:
    plt.plot(EFs, Ef_vs_EF[q], label=f"q={q}")
plt.xlabel("E_F (eV, relative to VBM)")
plt.ylabel("Formation energy (eV)")
plt.legend()
plt.show()
