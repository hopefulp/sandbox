#!/home/joonho/anaconda3/bin/python
import aselite
import re
import argparse
import os
def conv_pos2xyz(fname):
    ### read atoms
    if os.path.isdir(fname):
        atoms = aselite.read_vasp(fname+'/CONTCAR')
    else:
        atoms = aselite.read_vasp(fname)
    if re.match('POSCAR.', fname) or re.match('CONTCAR.', fname):
        fnamelist=fname.split('.')
        outf=fnamelist[1]+'.xyz'
    elif fname.endswith('POSCAR') or fname.endswith('CONTCAR'):
        if '/' in fname:
            fnlist=fname.split('/')
            outf=fnlist[-2]+'.xyz'
        else:
            outf=fname+'.xyz'
    ### if dirname
    elif os.path.isdir(fname):
        outf  = fname + '.xyz'
    print(f"output fname = {outf}")
    aselite.write_xyz('%s' % outf, atoms)
    return 0

def main():

    parser = argparse.ArgumentParser(description="change POSCAR|CONTCAR to xyz format  ")
    parser.add_argument('posname', help="pos2xyz.py [dir/CONTCAR|CONTCAR.name]")
    args = parser.parse_args()

    conv_pos2xyz(args.posname)

if __name__ == "__main__":
    main()



