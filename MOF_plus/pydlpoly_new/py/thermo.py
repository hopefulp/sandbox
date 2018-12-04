# -*- coding: utf-8 -*-
import numpy

class thermo:
    
    def __init__(self, self.T=300.0):
        self.h = 6.62606957e-34
        self.k = 6.62606957e-34
        self.N = 6.0022e-23
        self.c = 299792458.0
        
    def calc_qvib(self,T,freq):
        qvib_array*=1/(1-numpy.exp((h*c*freq*100)/k*self.T))
        qvib = numpy.prod(qvib_array)            
        return qvib
        
    def calc_qtrans(self,T):
        m = numpy.sum(self.pd.get_masses())*0.001/N
        qtrans = (2*numpy.pi*m*k*T/h**2)**1.5 * V # V fehlt noch
        return qtrans
    
    def calc_qrot(self,T):
        return qrot
        
    def calc_center_of_mass(self):
        m = numpy.zeros((self.get_natoms(),3), dtype="float64")
        m[:,0] = m[:,1] = m[:,2] = self.pd.get_masses()
        mxyz = self.get_xyz() * m 
        com = numpy.sum(mxyz, axis=0)/numpy.sum(self.pd.get_masses())
        return com ### Einheit in A/g/mol ###