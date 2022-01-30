#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os
from mod_doscar import obtain_doscar_head, change_Bheadline
def modify_doscar(fname, ext):
    of_save  = fname+ext
    of_tmp = fname+'tmp'
    fout=open(of_tmp, 'w')

    Nprelude = 5
    with open(fname, 'r') as f:
        ### analyze DOSCAR
        natom, Emax, Emin, ngrid, Ef, headline, E2nd = obtain_doscar_head(fname)

        ### change headline: ngrid-1, Emin-delE
        new_ngrid = int(ngrid) - 1
        new_headline = change_Bheadline(headline, Emin[:6], E2nd, ngrid, str(new_ngrid))
        totline = Nprelude + (natom+1) * (int(ngrid)+1) 

        for i, line in enumerate(f):
            if i < Nprelude:
                pass
                ### : if not, just fout.write(line)
            ### change block head line: Emin, ngrid
            elif line == headline:
                line = new_headline
            else:
                ene = line.strip().split()[0]
                ### skip and write twice each first line in the energy block
                if float(ene) == float(Emin):
                    continue
            ### just write if not error line
            fout.write(line)
    fout.close()
    ### change filename: DOSCAR --> DOSCAR_o, DOSCARtmp --> DOSCAR
    s1 = f"mv {fname} {of_save}"
    s2 = f"mv {of_tmp} {fname}"
    print(s1)
    print(s2)
    os.system(s1)
    os.system(s2)
    return 0

def main():
    parser = argparse.ArgumentParser(description='To get rid of abnormal DOS at start energy')
    parser.add_argument('file', nargs='?', default='DOSCAR', help='read DOSCAR')
    parser.add_argument('-o', '--old_suffix', default='_o', help='save original DOSCAR to DOSCAR_o')
    args = parser.parse_args()
    
    modify_doscar(args.file, args.old_suffix) 

if __name__ == '__main__':
    main()
