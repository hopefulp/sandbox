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
    parser.add_argument('-z', '--magmom', default=0.0, type=float, help='initial magnitude of magnetism')
    parser.add_argument('-i', '--ini', default='start', choices=['start','chg','wav','cont'],help='start/cont(chg/wav)')
    parser.add_argument('-e', '--cutoff', default=450, type=int, help='cut off energy')
    parser.add_argument('-a', '--precision', default='high', choices=['high','accurate'], help='precision')
    parser.add_argument('-r', '--relax', default='atom', choices=['sp', 'atom', 'cell'], help='atom/cell relaxation')
    parser.add_argument('-t', '--dft', default='pe', choices=['lda','gga','pe','rp','re','re0','pbe0','b3lyp','hse06','mk','ml'], help='gga method')
    parser.add_argument('-p', '--postscf', choices=['dos','pchg','band'], help='post-scf calculation')
    parser.add_argument('-d', '--dispersion', default='d3', choices=['d2', 'd3'], help='dispersion  scheme')
    parser.add_argument('-o', '--soc', action='store_true', help='whether soc or not')
    parser.add_argument('-u', '--uterm', type=float, help='input U-term energy')
    parser.add_argument('-l', '--log', default=1, type=int, choices={0,1,2,3}, help='{0:No,1:WAVECAR,2:CHGCAR,3:CHG')
    parser.add_argument('-v', '--solvent', action='store_true', help='whether solvent')
    parser.add_argument('-k', '--dielectric', default=78.3, type=float, help='dielectric constant')
    parser.add_argument('-md', '--dynamics', choices=['nvt', 'npt'], help='molecular dynamics flavor')
    args = parser.parse_args()

    ### 1. obtain default vasp repository
    arg_dic = vars(args) 
    #print arg_dic
    
    myvasp.make_incar(arg_dic, args.rw, args.iofile)


if __name__ == '__main__':
    main()
