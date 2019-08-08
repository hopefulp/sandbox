### energy conversion 
ev2kj = 96.4869
C100kj2ev = 1.036
hr2ev = 27.2107
j2cal = 0.23900
hr2kj = hr2ev * ev2kj
hr2_100kj = hr2ev * 0.964869 

class Q_Chem_aimd:
    def __init__(self, xyz='View.xyz', force='NucForces', energies='EComponents'):
        self.xyz=xyz
        self.force=force
        self.energies=energies

### Q-Chem model
Q_Chem_Epot = {'h2o': {'631gs': -76.4089462650}}

class Q_Chem_pot:
    def __init__(self, basis="631gs", mol='h2o'):
        self.basis=basis
        self.mol=mol
        self.pot=Q_Chem_Epot[self.mol][self.basis]
    # make this class callable
    def __call__(self):
        return self.pot
