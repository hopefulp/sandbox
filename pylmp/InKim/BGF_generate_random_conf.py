import sys
import getopt
import math
import random

import tqdm
import numpy as np

import bgf
import bgftools as bt
import scipy.spatial
import mathtools

version = '161115'

def query_safe_coords(exist, new, r=1.0):
    if not exist:
        return True # empty land

    T = scipy.spatial.KDTree(exist)
    for i in new:
        if T.query_ball_point(i, r=r):
            return False

    return True


def main(bgf_file, out_file, ff_file="", n=0, density=0.9, r=1.0):
    n_total_trial = 0
    # if molecule numbers are not specified, create the box with the existing molecules in the BGF file.
    # if molecule numbers specified, copy the molecules in the box as many as the number specified.
    if n:
        import copy
        mybgf = bgf.BgfFile()
        id = 0
        for i in tqdm.tqdm(range(int(n)), desc="Cloaning molecules", ncols=120):
            id += 1
            mybgf2 = bgf.BgfFile(bgf_file)
            for i in mybgf2.a:
                i.rNo = id
            mybgf = mybgf.merge(mybgf2)
    else:
        mybgf = bgf.BgfFile(bgf_file)

    # calculate cubic size
    mass = bt.getMass(mybgf, ff_file=ff_file)
    target_density = density;
    volume = mass / target_density / 6.022 * 10
    cell_x = math.pow(volume, 1.0/3)
    print("The script will generate a cubic with dimensions %f^3" % cell_x)
    mybgf.CRYSTX = [cell_x, cell_x, cell_x, 90.0, 90.0, 90.0]

    # get com
    redefined_coords = []
    molecules = bt.getMoleculeList(mybgf)
    for mol in tqdm.tqdm(molecules, ncols=120):
    #for index, mol in enumerate(molecules):
        n_trial = 0
        while True:
            n_trial += 1
            n_total_trial += 1
            cx, cy, cz = bt.getCom(mybgf, ff_file=ff_file, aNo_list=mol)
            new_pos = [random.uniform(0, cell_x) for i in range(3)]
            new_rot = [random.uniform(0, 2*math.pi) for i in range(3)]
            U = mathtools.rotate_matrix(new_rot[0], new_rot[1], new_rot[2])

            # check the room before lying
            mol_coords = []
            for ano in sorted(mol):
                atom = mybgf.getAtom(ano); x = atom.x; y = atom.y; z = atom.z
                x -= cx; y -= cy; z -= cz

                # rotate
                v = np.matrix([x, y, z]).T
                Uv = U * v
                x = float(Uv[0]); y = float(Uv[1]); z = float(Uv[2])

                # move CM to new position
                x += new_pos[0]; y += new_pos[1]; z += new_pos[2]

                mol_coords.append([x, y, z])

            # only successful trials can quit the while loop
            if query_safe_coords(redefined_coords, mol_coords, r=r):
                for index, ano in enumerate(sorted(mol)):
                    atom = mybgf.getAtom(ano)
                    atom.x = mol_coords[index][0]
                    atom.y = mol_coords[index][1]
                    atom.z = mol_coords[index][2]

                redefined_coords += mol_coords
                break
            else:
                # fail if too many trials performed
                if n_trial > 1000:
                    print("Failed to find a suitable coordinates to insert a molecule %s" % mol)
                    sys.exit(0)

    # check
    print("Checking bad contacts..")
    t = scipy.spatial.KDTree(redefined_coords)
    d = t.query_pairs(r)
    if not d:
        print("There are no bad contacts with distance < %.2f" % r)
    else:
        print("Looks there're %d bad contacts" % len(d))

    # save
    mybgf.saveBGF(out_file)
    print("File saved to %s" % out_file)
    print("\n** Stats **")
    print("Total generated coordinates: %d" % n_total_trial)
    print("Number of last trials: %d" % n_trial)
    print("Done.")


if __name__ == '__main__':

    usage = '''%s -b bgf_file -o out_file -f ff_file (-n n_molecules)''' % sys.argv[0]
    if len(sys.argv) < 2:
        print(usage);
        sys.exit(0)

    n_mol = 0
    radius = 1.0

    options, args = getopt.getopt(sys.argv[1:], 'hb:o:f:n:d:r:', ['help','bgf=','out=', 'nmolecules=', 'density=', 'radius='])
    for option, value in options:
        if option in ('-h', '--help'):
            print(usage)
            sys.exit(0);
        elif option in ('-b', '--bgf'):
            bgf_file = value
        elif option in ('-o', '--out'):
            out_file = value
        elif option in ('-f', '--ff'):
            ff_file = value
        elif option in ('-n', '--nmolecules'):
            n_mol = int(value)
        elif option in ('-d', '--density'):
            density = float(value)
        elif option in ('-r', '--radius'):
            radius = float(value)
        elif option == NULL:
            print(usage)
            sys.exit(0)

    # default
    if out_file == "":
        out_file = os.path.basename(bgf_file).split(".bgf")[0] + "_mod.bgf"

    # main call
    print(sys.argv[0] + " version " + str(version))
    main(bgf_file, out_file, ff_file=ff_file, n=n_mol, density=density, r=radius)

