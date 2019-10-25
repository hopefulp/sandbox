class data():
    def __init__(self, filename):
        assert isinstance(filename, str), "Filename for LAMMPS data is invalid."
        self.remark = ""
        self.n_atoms = 0
        self.n_bonds = 0
        self.n_angles = 0
        self.n_diheds = 0
        self.n_impropers = 0
        self.n_atom_types = 0
        self.n_bond_types = 0
        self.n_angle_types = 0
        self.n_dihed_types = 0
        self.n_improper_types = 0
        self.atoms = dict()
        self.bonds = dict()
        self.angles = dict()
        self.diheds = dict()
        self.bond_types = dict()
        self.angle_types = dict()
        self.dihed_types = dict()
        self.pair_types = dict()

        with open(filename, 'r') as f:
            lines = f.readlines()

        self.remark = lines[0]
