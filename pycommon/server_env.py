#### server name

import os

PPNs={'chi':1, 'login':4}
PBS_QChem_scripts={'chi':'None', 'login':''}
PBS_VASP_scripts={'chi':'None', 'login':''}

class Machines:
    """ info for KAIST server """
    def __init__(self, ppn, qc_pbs, vasp_pbs):
        self.ppn = ppn
        self.qc_pbs = qc_pbs
        self.vasp_pbs = vasp_pbs

_HOST = os.getenv('HOST')
_USER = os.getenv('USER')

if _HOST == 'chi':
    home = "/home/"+ _USER
elif _HOST == 'login':
    home = "/gpfs/home/"+ _USER
else:
    print("error in hostname %s" % _HOST)
python      = home + '/anaconda3/bin/python'
PBS_HOME    = home + '/sandbox/pypbs/'    

pbs_vname   = PBS_VASP_scripts[_HOST]
pbs_qname   = PBS_QChem_scripts[_HOST]
host_machine=Machines(PPNs[_HOST], PBS_HOME + pbs_qname, PBS_HOME + pbs_vname)



