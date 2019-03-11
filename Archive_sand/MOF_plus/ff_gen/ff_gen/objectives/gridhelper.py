import numpy
import string


class grid:
    """ generic grid generator class"""
    
    def __init__(self, mol, res = 0.2, offset = 4.0, scale = 1.0, sym = True):
        self.mol = mol
        self.res = res
        self.offset = offset
        self.scale = scale
        self.ngrid = 0
        self.grid = None
        self.grid_type = "generic"
        if sym == True:
            self.mol.symm.center()
            self.mol.symm.detect_symmetry()
        print self.mol.symm.symm
        return

    
    def within_vdw_sphere(self, point, scale=1.0):
        vdwr = self.mol.vdwr*scale
        r = self.mol.xyz-point
        d = numpy.sqrt(numpy.sum(r*r, axis=1))-vdwr
        # find if there are negative numbers
        nd = numpy.where(d<0)[0]
        if len(nd) == 0:
            return False
        return True
        
        
    def gen_grid(self):
        gxyz=[]
        #low = self.mol.xyz-self.mol.vdwr-self.offset
        low = self.mol.xyz-self.offset
        lowind = numpy.argmin(low,axis=0)
        lowx = low[lowind[0],0]
        lowy = low[lowind[1],1]
        lowz = low[lowind[2],2]
        low = numpy.array([lowx, lowy, lowz])
        #high = self.mol.xyz+self.mol.vdwr+self.offset
        high = self.mol.xyz+self.offset
        highind = numpy.argmax(high,axis=0)
        highx = high[highind[0],0]
        highy = high[highind[1],1]
        highz = high[highind[2],2]
        high = numpy.array([highx, highy, highz])
        self.li = numpy.array(numpy.fix(low/self.res)-1,"i")
        self.hi = numpy.array(numpy.fix(high/self.res)+2,"i")
        if self.mol.symm.symm.count("xy") : self.li[2]=0
        if self.mol.symm.symm.count("xz") : self.li[1]=0
        if self.mol.symm.symm.count("yz") : self.li[0]=0
        print self.li
        print self.hi
        for x in xrange(self.li[0],self.hi[0]):
            for y in xrange(self.li[1],self.hi[1]):
                for z in xrange(self.li[2],self.hi[2]):
                    point =numpy.array([x,y,z],"d")*self.res
                    gxyz.append(point.tolist())
                    if not self.within_vdw_sphere(point, scale=self.scale):
                        # yes .. keep this point
                        gxyz.append(point.tolist())
        self.grid = numpy.array(gxyz)
        self.ngrid = len(self.grid)
        return

class gridplot:

    def __init__(self, ngrid, grid, energies):
        assert numpy.shape(grid)[0] == ngrid
        assert ngrid = numpy.shape(energies)
        self.ngrid = ngrid
        self.grid = grid
        self.energies = energies
        return


    def get_xy(self):
        return

    def get_xz(self):
        return

    def get_yz(self):
        return


            
        
    
