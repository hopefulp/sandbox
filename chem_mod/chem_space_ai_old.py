""" Atom array for molecule """
ethylene    = ['C', 'C', 'H', 'H', 'H', 'H']
diss_cho    = ['C', 'O', 'H', 'H', 'H']
h2o         = [ 'O', 'H', 'H']
mol2atoms={'ethylene':ethylene, 'Diss_CHO':diss_cho, 'h2o':h2o}

atomic_number = { 'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10}

lattice_param_default = 10.0
lattice_param_h2o = (15.6404, 15.6404, 15.6404)
lattice = 'Lattice="{a:4.1f} 0.0 0.0 0.0 {a:3.1f} 0.0 0.0 0.0 {a:3.1f}" '.format(a=lattice_param_default)
lattice_ortho = 'Lattice=" {a:8.4.1f} 0.0 0.0 0.0 {a:8.4f} 0.0 0.0 0.0 {a:8.4f}" '
properties_e = 'Properties="species:S:1:pos:R:3:Z:I:1" '             # only position
properties_ef = 'Properties="species:S:1:pos:R:3:Z:I:1:forces:R:3" '   # include forces
lattice_properties_e=lattice+properties_e
lattice_properties_ef=lattice+properties_ef


class Bulk:
    def __init__(self, lattice_param):
        self.lattice_param=lattice_param
        #self.lattice_str = 'Lattice=" {:8.4f} 0.0 0.0 0.0 {:8.4f} 0.0 0.0 0.0 {:8.4f}" '.format(lattice_param[0],lattice_param[1],lattice_param[2])
        self.lattice_str = self.latt_ortho(self.lattice_param)

    def latt_ortho(self, latt):
        st =  'Lattice="'
        st += '{:8.4f} 0.0 0.0 '.format(latt[0])
        st += '0.0 {:8.4f} 0.0 '.format(latt[1])
        st += '0.0 0.0 {:8.4f}" '.format(latt[2])
        return st

class H2O(Bulk):
    def __init__(self,lattice_param=None, property_type=None):
        if lattice_param == None:
            lattice_param = lattice_param_h2o
        Bulk.__init__(self, lattice_param)
        if property_type == None:
            property_type = 2
        self.property_str=self.get_properties(property_type)
        self.extxyz_str=self.lattice_str + self.property_str

    def get_properties(self, pro_type):
        if pro_type == 2:
            proper=properties_ef
        return proper
        

