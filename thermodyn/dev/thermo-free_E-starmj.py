#!/usr/bin/python3

# Update: 2021 / 05 / 14
# Min Jong Noh (starnmj@kaist.ac.kr)

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

############    Functions    #########################################

def file_read(filename):
    lineinfo = []; wordinfo = []
    
    with open(filename) as f:
        for i, l in enumerate(f):
            line = l; word = line.split()
            lineinfo.append(line); wordinfo.append(word)

    return lineinfo, wordinfo

def get_info(filepath):
    """ 
    This function will automatically search
      1) VASP total energy from OUTCAR 
      2) Atom(s) information from CONTCAR (ver.5 format)
    """
    current_path = os.popen('pwd').read()
    path = filepath
    
    E = []; atom = []
  
    os.chdir(path)
    if os.path.isfile('OUTCAR'):
        lineinfo, wordinfo = file_read('OUTCAR')
        for i in range(len(lineinfo)):
            if 'y  w' in lineinfo[i]:
                TE = float(wordinfo[i][-1])
                E.append(TE)
            else:
                pass
    else:
        print("No OUTCAR in", path)
    
    # (Caution) CONTCAR file must be VASP ver.5 format
    # line 6th : atomic element
    # line 7th : # of atoms for each element
    
    if os.path.isfile('CONTCAR'):
        lineinfo, wordinfo = file_read('CONTCAR')
        atom_element = wordinfo[5]; atom_number = wordinfo[6]
        for j in range(len(atom_element)):
            info = {atom_element[j] : int(atom_number[j])}
            atom.append(info)
    else:
        print("No CONTCAR in", path)
    
    os.chdir(current_path.strip())
    optE = float(E[-1])
    
    return optE, atom

def format_data(legend, paths):
    data = []
    if len(legend) == len(paths):
        for i in range(len(legend)):
            optE, atom = get_info(paths[i])
            data.append([legend[i], optE, atom])
    else:
        print("number of legend and path are not matched")

    return data

############    Example: lowest O-termination in MoS2 nanoribbon #####


# Step 1: set the VASP output paths after geometry optimization calculation

""" 
You must modify this variables to calculate your own simulation
"""

legend_target = ["O-000", "O-050", "O-100", "S-100"]
legend_refer  = ["a-MoO3", "bcc-Mo", "alpha-S", "mol-O2"]

path_target = [
             "/home2/starnmj/TASK/35.MoS2-oxidation/thermodynamic-tutorial/01.target_MoS2_nanoribbon/01.O-000",
             "/home2/starnmj/TASK/35.MoS2-oxidation/thermodynamic-tutorial/01.target_MoS2_nanoribbon/02.O-050",
             "/home2/starnmj/TASK/35.MoS2-oxidation/thermodynamic-tutorial/01.target_MoS2_nanoribbon/03.O-100",
             "/home2/starnmj/TASK/35.MoS2-oxidation/thermodynamic-tutorial/01.target_MoS2_nanoribbon/04.S-100"
              ]

path_refer = [
             "/home2/starnmj/TASK/35.MoS2-oxidation/thermodynamic-tutorial/02.reference/bulk/01.alpha-MoO3",
             "/home2/starnmj/TASK/35.MoS2-oxidation/thermodynamic-tutorial/02.reference/bulk/02.bcc-Mo",
             "/home2/starnmj/TASK/35.MoS2-oxidation/thermodynamic-tutorial/02.reference/bulk/03.alpha-S",
             "/home2/starnmj/TASK/35.MoS2-oxidation/thermodynamic-tutorial/02.reference/mol/01.mol-O2",
             ]

data_target = format_data(legend_target, path_target)
data_refer = format_data(legend_refer, path_refer) 

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
    print(legend_target[i], rich_point[i], poor_point[i])

# Step 4: visualization with matplotlib

X_shift_max = mu_O_rich - 0.5*E_molO2
X_shift_min = mu_O_poor - 0.5*E_molO2

plt.figure(figsize=(7,5))
for i in range(len(legend_target)):
    line = plt.plot([X_shift_max, X_shift_min], [rich_point[i], poor_point[i]], label='%s' % legend_target[i])

plt.axis([X_shift_min, X_shift_max, min(rich_point), max(poor_point)])
plt.vlines(X_shift_min+0.5, min(rich_point), max(poor_point), color='red', linestyle=':')
plt.vlines(X_shift_max-0.5, min(rich_point), max(poor_point), color='red', linestyle=':')
plt.xlabel('O chemical potential [eV]', fontsize=10)
plt.ylabel('E$_{edge}$ [eV]', fontsize=10)
plt.legend()
plt.show()

################################################################
