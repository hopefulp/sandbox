### analysis of DOSCAR
from common import whereami
import sys

nheadline = 5
dos_err = 100

def obtain_doscar_head(fname):
    '''
    read DOSCAR
    return  int:    natom, ngrid
            float:  Emax, Emin, Ef
            str:    bheadline (block headline)
            boolean: check doserr exists
    '''
    with open(fname, 'r') as f:
        ### analyze DOSCAR
        for i, line in enumerate(f):
            if i < nheadline:
                ### natom natom 1 0
                if i == 0:
                    natom = int(line.strip().split()[0])
                ### : if not just fout.write(line)
            elif i == nheadline:
                bheadline = line
                info = line.strip().split()
                Emax    = float(info[0])
                Emin    = float(info[1])
                ngrid   = int(info[2])
                Ef      = float(info[3])
            ### check whether DOS at first energy has error
            elif i == nheadline +1:
                elist   = line.strip().split()
                dos1    = float(elist[1])
                ncol_tdos = len(elist)
            elif i == nheadline + 2:
                dos2    = float(line.strip().split()[1])
                break
    if dos2 * dos_err < dos1:
        Ldos_err = True
    else:
        Ldos_err = False
    if ncol_tdos == 3:
        Lspin = False
    elif ncol_tdos == 5:
        Lspin = True
    else:
        print(f"DOSCAR parsing error: ncol of tdos {ncol_tdos}")
        sys.exit(100)
    print(f"{whereami():>15}(): DOS err {Ldos_err} at 1st energy {dos1} then {dos2} in {__name__}.py")
    return natom, Emax, Emin, ngrid, Ef, bheadline, Ldos_err, Lspin

def change_Bheadline(old, Emin, Erep, ngrid, new_ngrid):
    print(f"{old}")
    new_grid = old.replace(ngrid, new_ngrid)
    print(f"{new_grid}")
    new_headline = new_grid.replace(Emin, Erep)
    print(f"{new_headline} with {Erep} for {Emin} ")
    return new_headline  
