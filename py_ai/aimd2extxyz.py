#!/home/joonho/anaconda3/bin/python

import argparse
import os
import sys
import subprocess
from my_chem import *
import chem_space as cs

#import CenterShift
#import mdene
#import mdForce

fqchems = Q_Chem_aimd()
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

def aimd2extxyz(job, fxyz, fene, ffor, Lattice_Size):

    outf = job + '.extxyz'

    atomic_number = cs.atomic_number
    ### Advance by natom+2
    lxyz_coord = fxyz.readlines()        # natom+2 lines per frame
    natom = int(lxyz_coord[0])
    nframe1 = cal_nframes(fqchems.xyz, 0, natom+2)
    print(f"{nframe1} frames in {fqchems.xyz}")
    ### Advance by 1
    lenergies = fene.readlines()
    nframe2 = cal_nframes(fqchems.energies, 1, 1) 
    print(f"{nframe2} frames in {fqchems.energies}")
    ### Advance by 1
    lforces = ffor.readlines()
    nframe3 = cal_nframes(fqchems.force, 1, 1)
    print(f"{nframe3} frames in {fqchems.force}")
    nframe = min(nframe1, nframe2, nframe3)


    i=0
    iframe=0
    with open(outf, 'w') as f:
        ### scan xyz file and modify it
        for xyzline in lxyz_coord:
            if i%(natom+2)==0:
                iframe+=1
                if iframe > nframe:
                    break
                iatom=0
                f.write(xyzline)
                i+=1
                continue
            elif i%(natom+2)==1:
                f.write('Lattice="{0:<.1f} 0.0 0.0 0.0 {0:<.1f} 0.0 0.0 0.0 {0:<.1f}" '.format(Lattice_Size))
                f.write('Properties="species:S:1:pos:R:3:Z:I:1:forces:R:3" ')
                E_pot_hr = float(lenergies[iframe].strip().split()[2])
                if iframe==1:
                    E_pot_hr0 = E_pot_hr
                    epot = 0.0  # in unit of 100 kJ/mol
                else:
                    epot = (float(E_pot_hr) - E_pot_hr0)*hr2kj
                f.write('energy={:<.10f}\n'.format(epot))
                i+=1
            else:
                f.write(xyzline.strip())
                list_atomline=xyzline.split()
                f.write('{0:>2d}   '.format(atomic_number[list_atomline[0]]))
                fx, fy, fz = obtain_Force(lforces[iframe], iatom)
                f.write('{:.8f} {:.8f} {:.8f} \n'.format(float(fx), float(fy), float(fz)))
                iatom+=1
                i+=1
                
            #for atom in range(numstep*numatoms,numstep*numatoms+numatoms):
            #    f.write('{0:<5s}{1:<15.10f}{2:<15.10f}{3:<15.10f}'.format(Atoms[atom],X[atom],Y[atom],Z[atom]))
            #    f.write('{0:>2d}   '.format(atomic_number[Atoms[atom]]))
            #    for direc in range(atom*3, atom*3+3):
            #        f.write('{:.8f}   '.format(ForceList[direc]))
            #    f.write('\n')
    f.close()
    return 0

def hack_qchem_AIMDdir(job, dir1, lattice_size):
    pwd =  os.getcwd()
    if not dir1:
        dir1 = pwd
    if os.path.isdir(dir1):
        if fqchems.xyz in os.listdir(dir1):
            dir_aimd = dir1
        elif 'AIMD' in os.listdir(dir1):
            dir_aimd = dir1+'/AIMD'
            #os.chdir(dir_aimd)
        else:
            print("there is no Q-Chem AIMD outfiles")
            sys.exit(1)
    fcoord  = dir_aimd + '/' + fqchems.xyz
    fenergy = dir_aimd + '/' + fqchems.energies
    fforce  = dir_aimd + '/' + fqchems.force

    fxyz = open(fcoord,'r')
    fene = open(fenergy, 'r')
    ffor = open(fforce, 'r')
    #Center_Coord = CenterShift.FindCenter(XYZfin, numatoms) 
    #Lattice_Size = CenterShift.BoxSize(XYZfin, numatoms, Center_Coord)
    #Atoms, X, Y, Z = CenterShift.NewCoord(XYZfin, numatoms, Center_Coord, int(args.steps))
    #EnergyList = mdene.ExEne(Energy,int(args.steps))
    #ForceList = mdForce.ExForce(NucFor,numatoms,int(args.steps))
   
    #aimd2extxyz(outfin, numatoms, Lattice_Size, EnergyList, Atoms, X, Y, Z, atomic_number, ForceList, int(args.steps)) 
    aimd2extxyz(job, fxyz, fene, ffor, lattice_size)
    fxyz.close()
    fene.close()
    ffor.close()
    return 0

def main():
    parser = argparse.ArgumentParser(description='make extxyz from qchem output file: AIMD')
    parser.add_argument('-f', '--fjob', default='test', help='job or filename for output')
    parser.add_argument('-d', '--onedir', help='directory where AIMD output files exist')
    parser.add_argument('-ls', '--lattice_size', default=10.0, type=float, help='size of lattice vector in cubic')
    #parser.add_argument('-n','--steps', help='the number of MD step')

    args = parser.parse_args()    

    hack_qchem_AIMDdir(args.fjob, args.onedir, args.lattice_size)

if __name__ == '__main__':
    main()
