#!/usr/bin/env python
# -*- coding: utf-8 -*-
  
def read_restart(self, filename):
    """
    read a complete population from a restart file 
    
    :Parameters:
        - filename : filename of the restart file (V1.0)
    """
    rf = open(filename, "r")
    l = string.split(rf.readline())
    assert l[4] == "V1.0" , "Wrong version of restart file"
    l = string.split(rf.readline())
    l = string.split(rf.readline())
    self.n   = int(l[1])
    self.npop  = int(l[3])
    ngen  = int(l[5])
     assert self.n  == n
     assert self.np == npop
    self.ngen = ngen
    fit = np.zeros([self.npop], dtype="float64")
    val = np.zeros([self.npop, self.n], dtype="float64")
    for i in xrange(self.npop):
        l = string.split(rf.readline())
        assert i+1 == int(l[0])
        fit[i] = float(l[2])
        l = string.split(rf.readline())
        v = np.array(map(string.atof, l))
        val[i,:] = v
    rf.close()
    self.val = val
    return


def analyze(self, n, ric, atoms):
    for i in range(n):
        params = self.val[i,:]
        self.pd.write_variables(params)
        self.obj()
        j, j_glob = self.obj.map_ric(ric, atoms)
        delt_force = (self.obj.force - self.obj.ref_force)*self.obj.fact_force
        both = np.zeros((np.shape(delt_force)[0],2))
        if ric == 'tor':
            sort = np.argsort(self.obj.d_dihedrals[:,j])
            for l in range(np.shape(sort)[0]):
                index = sort[l]
                both[l,0] = self.obj.d_dihedrals[index,j]
                both[l,1] = delt_force[index, j_glob]
        if ric == 'ibe':
            sort = np.argsort(self.obj.d_angles[:,j])
            for l in range(np.shape(sort)[0]):
                index = sort[l]
                both[l,0] = self.obj.d_angles[index,j]
                both[l,1] = delt_force[index, j_glob]
        if ric == 'str':
            sort = np.argsort(self.obj.d_bonds[:,j])
            for l in range(np.shape(sort)[0]):
                index = sort[l]
                both[l,0] = self.obj.d_bonds[index,j]
                both[l,1] = delt_force[index, j_glob]
        f = open('%s_%s_%s.dat' % (i, ric, atoms), 'w')
        for k in range(np.shape(delt_force)[0]):
            if ric == 'tor':
                f.write('%s %6.6f %6.6f\n' % (k, both[k,0], both[k,1]))
            if ric == 'ibe':
                f.write('%s %6.6f %6.6f\n' % (k, both[k,0], both[k,1]))
            if ric == 'str':
                f.write('%s %6.6f %6.6f\n' % (k, both[k,0], both[k,1]))
            if ric == 'obe':
                f.write('%s %6.6f\n' % (k, delt_force[k,j_glob]))
        f.close()
    return
 
