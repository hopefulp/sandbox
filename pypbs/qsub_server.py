#!/usr/bin/python

import argparse
import re
import os
from common import *

def print_sge(software):
    _HOME = os.getenv('HOME')
    if not software:
        print "Use:: qfree - to check freed node"
        print "      qhist - to check user occupancy"
        print "      qmem  - to check memory usage"
        print "      qstat - see my job"
        print "            -f for full"
        print "            -f | sed '/---/d' | sed 's*0/**' | sed 's*resv/**' "
        print "Usage:: -s software to see how to qsub qchem, grmx, etc"

    elif software == 'qchem':
        print "SGE(Rhee's group): Q-Chem"
        print "Usage:: qsub -N jobname -v job=qchem_input -v ver='3.2p' /gpfs/home/joonho/sandbox_gl/pypbs/sge_qchem.tcsh"
        print "        -N jobname: listed in 'qstat', jobname should be alphanumeric"
        print "        -v job: variable in 'qchem job.in job.out' in script sge_qchem.tcsh"
        print "        -v ver: version ['3.2p', '3.2s', '4.3s'] parallel and serial"
        q = "Do you want to see sge_qchem script?"
        if yes_or_no(q):
            sandbox = _HOME + '/sandbox_gl/pypbs/sge_qchem.tcsh'
            com = "more %s" % sandbox
            os.system(com)
    return 0

def job_description(server, software):
    if server == 'sge':
        print_sge(software) 


def main():
    parser = argparse.ArgumentParser(description='how to use qsub')
    parser.add_argument('server', default='sge', choices=['sge', 'chi'], help='jobname in pbs file')
    parser.add_argument('-s', '--software', choices=['qchem', 'grmx'], help='number of processor per node in the server')
    args = parser.parse_args()

    job_description(args.server, args.software) 

if __name__ == '__main__':
    main()
