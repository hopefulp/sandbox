#!/home/joonho/anaconda3/bin/python3

# Update: 2021/05/14  Min Jong  Noh     (starnmj@kaist.ac.kr)
# Update: 2022/08     Joonho    Park    (hopefulp@gmail.com)

"""
This example Python script provides post-processing calculation for constructing
edge (or Surface) free energy to find lowest-energy structure as a function of oxygen pressure.

It is only availabe for DFT calculations with VASP calculator.
Therefore, the users should write the individual VASP calculation file paths 
that sholud include output files as <OUTCAR, CONTCAR>.

See details in the referece below
- Composition, structure, and stability of RuO2(110) as a function of oxygen pressure (PRB 65, 035406, 2001)
"""

import os
import matplotlib.pyplot as plt
from thermo_env import data_target, data_refer, targets

############    Functions  -> rf. thermo_env   #########################################

############    Example: lowest O-termination in MoS2 nanoribbon #####


# Step 1: set the VASP output paths after geometry optimization calculation -> rf. thermo_env

# Step 2: get oxygen poor limit (0.5 delta_H) (cf. oxygen rich limit = 1/2E_O2)

""" 
Get enthalpy change of total reaction based on reference states
Insert exact stoichiometry of each reference state (fold_XXX)
"""

aMoO3 = data_refer[0] # a-MoO3
bccMo = data_refer[1] # bcc-Mo
molO2 = data_refer[3] # mol-O2
fold_aMoO3 = 4; E_aMoO3 = aMoO3[1]/fold_aMoO3
fold_bccMo = 2; E_bccMo = bccMo[1]/fold_bccMo
fold_molO2 = 1; E_molO2 = molO2[1]/fold_molO2

# Total reaction : delta_H = a-MoO3 - (1 * bcc-Mo + 3 * (1/2) * O2)
delta_H = E_aMoO3 - (1*E_bccMo + 3*0.5*E_molO2)

mu_O_rich = 0.5*E_molO2+0.5 # O rich limit
mu_O_poor = 0.5*E_molO2+0.5*delta_H-0.5 # O poor limit


# Step 3: get edge free energy (see PRB 67, 085410, 2003) & gradient

"""
edge free energy equation of O-termination MoS2 nanoribbon
         TE (system) - [n*mu(Mo)+m*mu(S)+l*mu(O)] (n, m, l: number of atoms in system)
G =   -----------------------------------------------------------------------------------
                 2L (2: both side, L: edge length)
"""
L = 3.16
aS = data_refer[2] # a-S
fold_aS = 128; E_aS = aS[1]/fold_aS

rich_point = []
poor_point = []
for i in range(len(data_target)):
    sys = data_target[i]; E_sys = sys[1]
    
    if len(sys[2]) == 3:
        n = sys[2][0]['Mo']
        m = sys[2][1]['S']
        l = sys[2][2]['O']
    elif len(sys[2]) == 2:
        n = sys[2][0]['Mo']
        m = sys[2][1]['S']
        l = 0
    
    G_rich = E_sys - (n*E_bccMo + m*E_aS + l*mu_O_rich)
    G_rich_norm = G_rich / (2 * L)
    rich_point.append(G_rich_norm)
    G_poor = E_sys - (n*E_bccMo + m*E_aS + l*mu_O_poor) 
    G_poor_norm = G_poor / (2 * L)
    poor_point.append(G_poor_norm)

for i in range(len(rich_point)):
    print(targets[i], rich_point[i], poor_point[i])

# Step 4: visualization with matplotlib

X_shift_max = mu_O_rich - 0.5*E_molO2
X_shift_min = mu_O_poor - 0.5*E_molO2

plt.figure(figsize=(7,5))
print(f"num legend target: {len(targets)}")
print(f"{targets}")
for i in range(len(targets)):
    line = plt.plot([X_shift_max, X_shift_min], [rich_point[i], poor_point[i]], label='%s' % targets[i])

plt.axis([X_shift_min, X_shift_max, min(rich_point), max(poor_point)])
plt.vlines(X_shift_min+0.5, min(rich_point), max(poor_point), color='red', linestyle=':')
plt.vlines(X_shift_max-0.5, min(rich_point), max(poor_point), color='red', linestyle=':')
plt.xlabel('O chemical potential [eV]', fontsize=10)
plt.ylabel('E$_{edge}$ [eV]', fontsize=10)
plt.legend()
plt.show()

################################################################
