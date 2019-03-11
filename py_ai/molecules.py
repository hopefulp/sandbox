""" Atom array for molecule """
ethylene = ['C', 'C', 'H', 'H', 'H', 'H']

mol2atoms={'ethylene':ethylene}

atomic_number = { 'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10}


lattice_param = 10.0
lattice = 'Lattice="{a:4.1f} 0.0 0.0 0.0 {a:3.1f} 0.0 0.0 0.0 {a:3.1f}" '.format(a=lattice_param)
properties_e = 'Properties="species:S:1:pos:R:3:Z:I:1" '             # only position
properties_ef = 'Properties="species:S:1:pos:R:3:Z:I:1:forces:R:3" '   # include forces
lattice_properties_e=lattice+properties_e
lattice_properties_ef=lattice+properties_ef
