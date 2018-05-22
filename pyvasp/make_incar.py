#!/usr/bin/python

import argparse
import myvasp
    


def main():
    
    global ini_dvasp, pwd
    parser = argparse.ArgumentParser(description='make incar')
    parser.add_argument('-f', '--iofile', default='incar.key', help='iofile for all arguments')
    parser.add_argument('-w', '--rw', default='w', help='read or write')
    parser.add_argument('-m', '--mag', default='nm', choices=['nm','fm','afm'], help='ON/OFF magnetism')
    parser.add_argument('-n', '--nmag', type = int, help='number of magnetic atoms')
    parser.add_argument('-z', '--magmom', default='3.0', help='initial magnitude of magnetism')
    parser.add_argument('-i', '--ini', default='start', choices=['start','chg','wav','cont'],help='start/cont(chg/wav)')
    parser.add_argument('-e', '--cutoff', default=450, type=int, help='cut off energy')
    parser.add_argument('-a', '--precision', default='high', choices=['high','accurate'], help='precision')
    parser.add_argument('-c', '--crelax', default='cell', choices=['sp', 'atom', 'cell'], help='atom/cell relaxation')
    parser.add_argument('-t', '--dft', default='pe', choices=['lda','gga','pe','pbe0','b3lyp','hse06','mk','ml'], help='gga method')
    parser.add_argument('-p', '--postscf', choices=['dos','pchg','band'], help='post-scf calculation')
    parser.add_argument('-d', '--dispersion', default='d2', choices=['d2', 'd3'], help='dispersion  scheme')
    parser.add_argument('-o', '--soc', action='store_true', help='whether soc or not')
    parser.add_argument('-u', '--uterm', type=float, help='input U-term energy')
    parser.add_argument('-l', '--log', default='True', help='|No|True|Full|, T writes wavecar & chgcar, full writes AECHG')
    parser.add_argument('-v', '--solvent', action='store_true', help='whether solvent')
    parser.add_argument('-k', '--dielectric', default='78.3', help='dielectric constant')
    args = parser.parse_args()

    ### 1. obtain default vasp repository
    arg_dic = vars(args) 
    #print arg_dic
    
    myvasp.make_incar(arg_dic, args.rw, args.iofile)


if __name__ == '__main__':
    main()
