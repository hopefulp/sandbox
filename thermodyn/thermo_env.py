import os


# Step 1: set the VASP output paths after geometry optimization calculation

""" 
You must modify this variables to calculate your own simulation
"""

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
      2) dict{atom species: natom} from CONTCAR
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
    '''
    return data
        list of legend, opt_energy, atom_format=dict{atom species:natoms}
    '''
    data = []
    if len(legend) == len(paths):
        for i in range(len(legend)):
            optE, atom = get_info(paths[i])
            data.append([legend[i], optE, atom])
    else:
        print("number of legend and path are not matched")

    return data


job_dir = '/home/joonho/research/MoS2oxi/VASP/thermodynamic-tutorial/DFT/'
dtarget = f"{job_dir}01.target_MoS2_nanoribbon/"
drefer  = f"{job_dir}02.reference/"

targets = ["O-000", "O-050", "O-100", "S-100"]
refers  = ["a-MoO3", "bcc-Mo", "alpha-S", "mol-O2"]

path_target = [	f"{dtarget}01.{targets[0]}", f"{dtarget}02.{targets[1]}",
               	f"{dtarget}03.{targets[2]}", f"{dtarget}04.{targets[3]}"
              ]

path_refer = [ 	f"{drefer}bulk/01.{refers[0]}", f"{drefer}bulk/02.{refers[1]}",
             	f"{drefer}bulk/03.{refers[2]}", f"{drefer}mol/01.mol-O2",
             ]

data_target = format_data(targets, path_target)
data_refer = format_data(refers, path_refer)

