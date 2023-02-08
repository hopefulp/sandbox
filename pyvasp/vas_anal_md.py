#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os
import subprocess

def anal_oszicar(fname, outf):
    if os.path.isdir(fname):
        f = fname + '/OSZICAR'
    else:
        f = fname
    st = f"grep T= {f} > mdlog.dat"
    subprocess.check_output(st, shell = True)

    t = []
    T = []
    Etot = []
    F = []
    Epot = []
    Ekin = []
    ther_pot = []
    ther_kin = []
    with open("mdlog.dat", 'r') as f:
        for line in f: 
            #print(line)
            li = line.strip().split()
            t.append(int(li[0]))
            T.append(int(float(li[2])))
            Etot.append(float(li[4]))
            F.append(float(li[6]))
            Epot.append(float(li[8]))
            Ekin.append(float(li[10]))
            ther_pot.append(float(li[12]))
            ther_kin.append(float(li[14]))
            print(f"{int(li[0]):>5} {int(float(li[2])):>5} {float(li[4]):12.5f} {float(li[6]):12.5f} {float(li[8]):12.5f} {float(li[10]):10.5f} {float(li[12]):10.5f} {float(li[14]):10.5f}")
            
    ### write to file
    fout=open(outf, 'w')
    for i in range(len(t)):
        fout.write(f"{t[i]:>5} {T[i]:>5} {Etot[i]:12.5f} {Epot[i]:12.5f} {Ekin[i]:10.5f}\n")
        ### check energy
        print(f"Etot {Etot[i]:12.5f} EF {F[i]:12.5f}  Epot+Ekin {Epot[i]+Ekin[i]:12.5f} Eall {Epot[i]+Ekin[i]+ther_pot[i]+ther_kin[i]:12.5f}")
    fout.close()

    return 0

def main():
    parser = argparse.ArgumentParser(description='read vasp logfile')
    parser.add_argument('file', help='read logfile or OSZICAR')
    parser.add_argument('-of', '--outf', default='md.dat', help='save md data')
    args = parser.parse_args()
    
    anal_oszicar(args.file, args.outf) 

if __name__ == '__main__':
    main()
