# Barometer which can apply a constant force into one direction on a molecule
# Allows system to adjust density during NVT simulation
# D. Pache (RUB 2015)

import numpy as np


class forcebar:
    
    def __init__(self, pd, forceconst, molecule, zref):
        self.pd = pd
        self.forceconst = forceconst
        self.molecule = self.pd.molecules.mlist[molecule]
        self.zref = zref
        self.force = np.zeros([self.molecule.natoms,3], "d")
        self.force[:,2] = -self.forceconst
        # for debug
        self.pd.pprint("Initalizing molecular barostat for molecule %d (%s)" % (molecule, self.molecule.mname))
        # self.out = open("forcebar.out", "w")
        return
    
    def __call__(self):
        # compute energy
        xyz = self.molecule.get_xyz()
        E   = np.sum(xyz[:,2]-self.zref)*self.forceconst
        self.molecule.add_force(self.force)
        #
        # self.out.write("Force:\n%s\n" % np.array2string(self.force[:10,2]))
        # self.out.write("XYZ:\n%s\n" % np.array2string(xyz[:10,2]))
        # self.out.flush()
        return E