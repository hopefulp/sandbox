#### to modify vasp input file
import sys
import vasp_job

def get_atoms_poscar(line):
    if atom in lines[0]:
        alist = lines[0].split()
    return alist

def fixedMD_POSCAR(poscar, atom, atoms=None):
    '''
    poscar  to be modified
    atom    kind to be moved for zpe calculation
    '''
    with open(poscar, 'r') as f:
        lines = f.readlines()
    with open('POSCAR', 'w') as f:
        iatom = 0
        for i, line in enumerate(lines):
            if i == 0:
                f.write(line)
                continue
            elif i < 5:
                f.write(line)
            ### atom lines could exist or not
            elif i==5 and any(j.isalpha() for j in line):
                atoms = line.strip().split()
                f.write(line)
            elif i <= 6 and any(j.isdigit() for j in line):
                natoms = list(map(int, line.strip().split()))
                if len(atoms) != len(natoms):
                    sys.exit(0)
                f.write(line)
                ### calculate
                ntotal = sum(natoms)
                ind = atoms.index(atom)
                nzpe = natoms[ind]
                isum = 0
                for i, na in enumerate(natoms):
                    if i < ind:
                        isum += na 
                f.write("Selective dynamics\n")
                print(f"ind {ind} pre sum {isum} ")
            elif i <= 7 and any(i.isalpha() for i in line):
                f.write(line)
            else:
                if iatom < isum:
                    new_line = line.rstrip() + " F F F\n"
                    f.write(new_line)
                elif iatom < isum + nzpe:
                    new_line = line.rstrip() + " T T T\n"
                    f.write(new_line)
                elif iatom < ntotal:
                    new_line = line.rstrip() + " F F F\n"
                    f.write(new_line)
                else:
                    break
                iatom += 1

        
        
    return


