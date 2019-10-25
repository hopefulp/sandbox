#!/home/joonho/anaconda3/bin/python
import argparse
from chem_space import *

### density calculation
def density2cube_by_n(mname, n, d):
    if not mname in MW.keys():
        print(f"add {mname} in chem_space")
        return 0
    else:
        mw = MW[mname]
        nmolecule_g=1./mw * NA          # 1./mw*NA; number in 1g
        vol_g = 1./d/ANG3               # 1./d/ANG3; vol of 1g in Ang^3
        vol   = vol_g/nmolecule_g*n     # vol_g/nmolecule_g*n; vol for 1 mol * n mol
        a_cubic = vol**(1./3.)
        print(f"cube for {n} molecule : {a_cubic:10.5f} {a_cubic:10.5f} {a_cubic:10.5f}")
        return 0

def main():
    parser = argparse.ArgumentParser(description="math for chem")
    parser.add_argument('-d','--density', type=float,  help="input density ")
    parser.add_argument('-n','--number', type=int, help="input density ")
    parser.add_argument('-m','--molecule', help="input molecule")
    args = parser.parse_args()

    density2cube_by_n(args.molecule, args.number, args.density)

if __name__ == "__main__":
    main()
