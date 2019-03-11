""" 
             lbfgs
             
     This is a python wrapper for the f90ified version of the original LBFGS code 
     by Jorge Nocedal
     http://users.eecs.northwestern.edu/~nocedal/lbfgs.html
     
     References:
     J. Nocedal. Updating Quasi-Newton Matrices with Limited Storage (1980), Mathematics of Computation 35, pp. 773-782.
     D.C. Liu and J. Nocedal. On the Limited Memory Method for Large Scale Optimization (1989), Mathematical Programming B, 45, 3, pp. 503-528.

     R. schmid RUB 2012
     
     Note: This is to avoid dependencies (avoid Scipy LBFGS)
     
     NOTE (2015): with adding MExt to allow multiple instances we need to pass the parent pd
     
"""
import numpy
from molsys import mpiobject
#import _pydlpoly
#dlp_lbfgs = _pydlpoly.lbfgs_module

    
class lbfgs(mpiobject):
    
    def __init__(self,pd,n,m,fc,thresh=0.1,verbose = True, mpi_comm = None, out = None):
        """ on init, the dimensions (N) must be specified, as well as the number of previous gradients to be 
            taken into account (M). The optimizer expects a flat numpy double array.
            upon instantiation the work array is allocated on f90 level.
            also a callback function providing the arrays needs to be given
            
            fc should be called with the flat input arrays x and g of length N, it should return the function value
            The trick is really that we sue the "fortran" style and pass on a numpy array as a "pointer" to the data.
            So the arrays x and g should be allocated on python level in the calling environment. after the function returns it
            contains the final function values.
            All computations are done "in place" to avoid extensive copying
            
            Note on IFLAG: this needs to be zero dim numpy array (not just a scalar integer)
            because its value is changed on fortran level.
        """
        super(lbfgs, self).__init__(mpi_comm,out)
        # NEW for MExt
        self.pd = pd
        self.dlp_lbfgs = self.pd._pydlpoly.lbfgs_module
        #
        self.n = n
        self.m = m
        self.fc = fc
        self.thresh = thresh
        self.dlp_lbfgs.lbfgs_init(self.n,self.m,self.thresh)
        self.iflag = numpy.array([0],"i")
        # output control
        self.write_every = 10
        self.verbose = verbose
        return

    def __del__(self):
        self.dlp_lbfgs.lbfgs_free()
        
    def __call__(self, x, g, maxiter= None, thresh=None):
        if thresh:
            self.thresh=thresh
        self.dlp_lbfgs.EPS = self.thresh
        f = self.fc(x, g)
        self.dlp_lbfgs.lbfgs(f,g,x,self.n,self.iflag)
        icount = 1
        stop = False
        while not stop:
            f = self.fc(x, g)
            self.dlp_lbfgs.lbfgs(f,g,x,self.n,self.iflag)
            if ((((icount-1)%self.write_every) == 0) and (self.verbose == True)):
                self.pprint("   LBFGS %4d: f= %12.6f rmsg= %12.6f" % (icount, f, self.dlp_lbfgs.RMSG))
            status = self.iflag[0]
            icount += 1
            if (status!=1) and (icount > 1): stop = True
            if maxiter:
                if icount > maxiter: stop = True
        self.pprint("   LBFGS %4d: f= %12.6f rmsg= %12.6f" % (icount, f, self.dlp_lbfgs.RMSG))
        if status == 0:
            self.pprint("   LBFGS is converged")
        else:
            # better error checking ... maybe reset and start again?
            self.pprint("   LBFGS not converged!!!!!!!")
            self.pprint(self.dlp_lbfgs.INFO)
        return
            
    def set_iprint(self, p1, p2):
        self.dlp_lbfgs.IPRINT[0] = p1
        self.dlp_lbfgs.IPRINT[1] = p2
        return

    def set_diag(self, diag):
        if len(diag) != self.n:
            raise ValueError, "Wrong length of diag array"
        diag = numpy.array(diag,"d")
        self.pprint(self.dlp_lbfgs.DIAGCO)
        self.dlp_lbfgs.DIAGCO = 1
        self.dlp_lbfgs.DIAG[:] = diag[:]
        self.pprint("using a diagonal hessian input")
        return
           


