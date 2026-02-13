#!/home/joonho/anaconda3/bin/python
'''
structout2fdf의 Docstring
read Siesta outfile and write STRUCT.fdf input file for next job
    Sesta outfile: Label.STRUCT_OUT
    Usage {file} Label.STRUCT_OUT [-o output_filename]
        output_filename will be STRUCT.fdf for Siesta input in the next job
'''
import argparse
import sys, os
from oorinano.calculator import siesta as s2

def write_siesta_struct(inf, outf):
    atom = s2.read_struct_out(inf)
    s2.writeAtomicStructure(atom, fname=outf)
    return 0

def main():
    parser = argparse.ArgumentParser(description='To get rid of abnormal DOS at start energy')
    parser.add_argument('file',  nargs='?', help='read atoms file in fdf')
    parser.add_argument('-o', '--outfile', default='struct_out.fdf', help='save original DOSCAR to DOSCAR_o')
    parser.add_argument('-u', '--usage', action='store_true')
    args = parser.parse_args()

    if args.usage:
        print(__doc__.format(file=os.path.basename(__file__)))
        sys.exit(0)
    write_siesta_struct(args.file, args.outfile)

if __name__ == '__main__':
    main()

