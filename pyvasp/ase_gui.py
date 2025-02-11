#!/home/joonho/anaconda3/bin/python
from ase import *
#from ase.calculators.vasp import Vasp


def vis(inp):
    atoms = inp
    view(atoms)


def main():
    parser = argparse.ArgumentParser(description='How to make a continuous job dir')
    parser.add_argument('-i', '--input', help='conceivable file by ase.io.read()')
    #parser.add_argument('-io', '--incar_option', nargs='*', help='input key-value pairs in the list from command line, if value=None, del')
    
    args = parser.parse_args()

    vis(args.input)
    return 0



if __name__ == 'main':
    main()
