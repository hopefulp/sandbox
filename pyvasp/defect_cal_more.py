#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# ---------------------
# 1) 입력: (예시) DFT 총에너지 값들 (eV)
# ---------------------
E_bulk = -2000.0               # pristine supercell total energy (eV)
# defect total energies for q=0,+1,+2 (toy values)
E_def = {0: -1996.0, 1: -1994.0, 2: -1991.0}

# chemical potentials (O2 reference)
E_O2 = -9.86                   # DFT energy of O2 (example)
muO_ref = 0.5 * E_O2           # reference mu_O
muO_Orich = muO_ref            # O-rich
muO_Op = muO_ref - 2.0         # O-poor (example: 2 eV lower)

# choose chemical condition
muO = muO_Orich                # change to muO_Op for O-poor

# defect stoichiometry: oxygen removed => n_O = -1
n_O = -1

# band parameters
Ev = 0.0                       # VBM reference
Eg = 3.0                       # bulk band gap (eV)
EFs = np.linspace(0, Eg, 201)

# ---------------------
# 2) Uncorrected formation energies
# Ef_uncorr = E_def - E_bulk - sum n_i mu_i + q*(E_F + Ev)
# ---------------------
Ef_uncorr = {}
for q, E_d in E_def.items():
    Ef0 = E_d - E_bulk - n_O * muO + q * (0 + Ev)   # formation energy at E_F=0
    # actual E_f(EF) = Ef0 + q * E_F
    Ef_uncorr[q] = Ef0

# ---------------------
# 3) Makov-Payne correction (leading term)
# Use MP: E_MP = q^2 * alpha * (e2) / (2 * eps * L)
# We'll implement a simple form with e2/(4*pi*eps0)=14.3996 eV·Å
# => E_MP(eV) = q^2 * alpha * 14.3996 / (2 * eps * L(Å))
# ---------------------
alpha = 2.837      # cubic Madelung constant (example)
eps_r = 100.0      # dielectric constant (example; TiO2 ~ high, pick 100 for demo) 
L = 20.0           # effective supercell length in Å (toy value)
conv = 14.3996     # e^2/(4*pi*eps0) in eV·Å

def E_mp(q):
    return (q*q) * alpha * conv / (2.0 * eps_r * L)

# ---------------------
# 4) FNV-like correction (illustrative)
# FNV is more involved; as a simple demonstration we take FNV_corr = f * MP_corr
# with f ~ 0.5-0.8 (FNV often smaller than naive MP for small cells)
# ---------------------
f_factor = 0.6
def E_fnv(q):
    return f_factor * E_mp(q)

# ---------------------
# 5) Build Ef(EF) arrays (uncorr, MP-corrected, FNV-corrected)
# ---------------------
Ef_uncorr_arr = {}
Ef_MP_arr = {}
Ef_FNV_arr = {}
for q in E_def.keys():
    Ef_uncorr_arr[q] = Ef_uncorr[q] + q * EFs
    Ef_MP_arr[q] = Ef_uncorr_arr[q] - E_mp(q)
    Ef_FNV_arr[q] = Ef_uncorr_arr[q] - E_fnv(q)

# ---------------------
# 6) Plot
# ---------------------
plt.figure(figsize=(7,5))
colors = {0:'C0', 1:'C1', 2:'C2'}
for q in Ef_uncorr_arr:
    plt.plot(EFs, Ef_uncorr_arr[q], label=f'q={q} (uncorr)', color=colors[q], linestyle='-')
    plt.plot(EFs, Ef_MP_arr[q], label=f'q={q} (MP corr)', color=colors[q], linestyle='--')
    plt.plot(EFs, Ef_FNV_arr[q], label=f'q={q} (FNV corr)', color=colors[q], linestyle=':')

# show transitions for one correction type (FNV)
# compute intersections pairwise for FNV-corrected lines
def intersection(E0_q, q1, E0_p, q2):
    # E_f = E0 + q*EF, intersection when E0_q + q1*EF = E0_p + q2*EF
    # EF = (E0_p - E0_q) / (q1 - q2)
    return (E0_p - E0_q) / (q1 - q2)

# get E0 at EF=0 for FNV-corrected
E0_fnv = {q: Ef_uncorr[q] - E_fnv(q) for q in E_def.keys()}
# intersections q=0 & q=1, q=1 & q=2
eps_01 = intersection(E0_fnv[0], 0, E0_fnv[1], 1)
eps_12 = intersection(E0_fnv[1], 1, E0_fnv[2], 2)

plt.axvline(eps_01, color='k', linestyle='--', alpha=0.6)
plt.text(eps_01+0.02, plt.ylim()[1]*0.7, f'ε(0/+1)={eps_01:.2f} eV', rotation=90)
plt.axvline(eps_12, color='k', linestyle='-.', alpha=0.6)
plt.text(eps_12+0.02, plt.ylim()[1]*0.5, f'ε(+1/+2)={eps_12:.2f} eV', rotation=90)

plt.xlim(0, Eg)
plt.ylim(-1, 6)
plt.xlabel('Fermi level E_F (eV, rel VBM)')
plt.ylabel('Formation energy E_f (eV)')
plt.title('Toy: V_O formation energy vs E_F (uncorr / MP / FNV)')
plt.legend(fontsize=8, ncol=2)
plt.grid(True)
plt.show()

# ---------------------
# 7) Print some numbers
# ---------------------
print("Parameters:")
print(f"  muO_ref = {muO_ref:.3f} eV (O-rich), muO chosen = {muO:.3f} eV")
print(f"  MP prefactor for q=1: {E_mp(1):.4f} eV, for q=2: {E_mp(2):.4f} eV")
print(f"  FNV (f={f_factor}) for q=1: {E_fnv(1):.4f} eV, for q=2: {E_fnv(2):.4f} eV")
print()
for q in sorted(E_def.keys()):
    print(f"q={q}: Ef(at EF=0) uncorrected = {Ef_uncorr[q]:.4f} eV, MP_corr={E_mp(q):.4f}, FNV_corr={E_fnv(q):.4f}")
print()
print(f"FNV transition levels: ε(0/+1)={eps_01:.3f} eV, ε(+1/+2)={eps_12:.3f} eV")
