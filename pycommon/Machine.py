#### server name
Servers = ('psi', 'kdft', 'mu', 'rho', 'chi')

PPNs={'kdft':12, 'psi':16, 'rho': 16, 'mu':32, 'chi':1}
PBS_QChem_scripts={'kdft': 'pbs_qchemkdft.sh', 'psi': 'pbs_qchem.sh', 'rho': 'pbs_qchem.sh', 'mu': 'pbs_qchemmu.sh'}

PBS_VASP_scripts={'kdft': 'pbs-vasp-kdft.sh', 'psi': 'pbs-vasp-psi.sh', 'rho': 'pbs-vasp-rho.sh', 'mu': 'pbs-vasp-mu.sh'}


class Machines:
    """ info for KAIST server """
    def __init__(self, ppn, qc_pbs, vasp_pbs):
        self.ppn = ppn
        self.qc_pbs = qc_pbs
        self.vasp_pbs = vasp_pbs

kdft=



