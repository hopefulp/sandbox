#!/home/joonho/anaconda3/bin/python
# modified from aimd2extxyz

import argparse
import my_chem
import chem_space as cs
from common import fname_decom
import sys
#import CenterShift
#import mdene
#import mdForce

def cal_nframes(fname, ntitle, n):
    com = f'wc -l {fname}'
    results = subprocess.check_output(com, shell=True)
    nline = int(results.split()[0])
    nframe = (nline - ntitle)/n
    return nframe


def obtain_Force(l_forces, iatom):
    list_force = l_forces.strip().split()
    ind = iatom*3 + 1
    return list_force[ind], list_force[ind+1], list_force[ind+2]

### lattice_size for cubic
def xyz2ext(f, lattice_size, energy, force=None):
    #aimd2extxyz(outfin, numatoms, Lattice_Size, EnergyList, Atoms, X, Y, Z, atomic_number, ForceList, int(args.steps)) 
    #def aimd2extxyz(job, dname, fxyz, fene, ffor, Lattice_Size):
    scale = 'hr2ev'
    escale = my_chem.__dict__[scale]

    f_pre, fext = fname_decom(f)
    if fext != "xyz":
        print(f"input file {f} does not have xyz extension")
        sys.exit(1)
    fp_xyz  = open(f, 'r')
    outf = f_pre + '.extxyz'

    atomic_number = cs.atomic_number
    ### Advance by natom+2
    lxyz_coord = fp_xyz.readlines()        # natom+2 lines per frame
    natom = int(lxyz_coord[0])
    #nframe1 = cal_nframes(fqchems.xyz, 0, natom+2)
    nframe = len(lxyz_coord)/(natom+2)
    #print(f"{nframe1} frames in {fqchems.xyz}")
    ### USE for many frames
    ### Advance by 1        
    #lenergies = fene.readlines()
    #nframe2 = cal_nframes(fqchems.energies, 1, 1) 
    #nframe2 = len(lenergies)-1
    #print(f"{nframe2} frames in {fqchems.energies}")
    ### Advance by 1
    #lforces = ffor.readlines()
    #nframe3 = cal_nframes(fqchems.force, 1, 1)
    #nframe3 = len(lforces)
    #print(f"{nframe3} frames in {fqchems.force}")
    #nframe = min(nframe1, nframe2, nframe3)

    i=0
    iframe=0
    #print(f"write to {outf}")
    with open(outf, 'w') as f:
        ### scan xyz file and modify it
        for xyzline in lxyz_coord:
            if i%(natom+2)==0:
                iframe+=1
                if iframe > nframe:
                    break
                iatom=0
                f.write(xyzline)
            elif i%(natom+2)==1:
                f.write('Lattice="{0:<.1f} 0.0 0.0 0.0 {0:<.1f} 0.0 0.0 0.0 {0:<.1f}" '.format(lattice_size))
                f.write('Properties="species:S:1:pos:R:3:Z:I:1:forces:R:3" ')
                #E_pot_hr = float(lenergies[iframe].strip().split()[2])
                epot = energy*escale

                f.write('energy={:<.10f}\n'.format(epot))
            else:
                f.write(xyzline.strip())
                list_atomline=xyzline.split()
                f.write('{0:>2d}   '.format(atomic_number[list_atomline[0]]))
                if not force:
                    fx, fy, fz = 0., 0., 0.
                else:
                    fx, fy, fz = obtain_Force(lforces[iframe], iatom)
                f.write('{:.8f} {:.8f} {:.8f} \n'.format(float(fx), float(fy), float(fz)))
                iatom+=1
            i+=1
                
    fp_xyz.close()
    print(f"output file is {outf} in Energy scale of {scale}")
    return 0

def main():
    parser = argparse.ArgumentParser(description='make extxyz from xyz format')
    parser.add_argument('fxyz', help='xyz file')
    parser.add_argument('-f', '--force', type=float, help='force file')
    parser.add_argument('-e', '--energy', type=float, help='force file')
    parser.add_argument('-ls', '--lattice_size', default=10.0, type=float, help='size of lattice vector in cubic')

    args = parser.parse_args()    

    xyz2ext(args.fxyz, args.lattice_size,args.energy, args.force)

if __name__ == '__main__':
    main()
