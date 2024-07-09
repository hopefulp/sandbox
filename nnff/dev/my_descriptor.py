### train descriptor
### Gaussian Descriptor

import numpy as np
from amp.descriptor.gaussian2 import make_symmetry_functions
import sys

class GS_param:
    '''
    func_param: 'log10', powNN=N^(m/N)
    Now eta & Rsm is calculated without Rc
    '''
    name='gs'
    def __init__(self, func_param='log10', pmin=0.05, pmax=5, pnum=4):
        self.func_param = func_param
        self.pmin = pmin
        self.pmax = pmax
        self.nparam = pnum

    def make_Gs(self, image):
        ### obtain atom symbols
        elements = set([atom.symbol for atom in image])
        elements = sorted(elements)     # because set is not subscriptable
        G = {}
        Np = self.nparam
        #etas = np.logspace(np.log(0.5), np.log(5), num=6)
        if self.func_param == 'log10':
            etas = np.logspace(np.log10(self.pmin), np.log10(self.pmax), num=Np)
        elif self.func_param == "powNN":
            ### G2a: for N, num(m) = N+1, m={0,1,...,N)
            etas = [ Np**(2*m/Np) for m in np.arange(Np+1)]
            Rs = list(np.zeros(len(etas)))
            ### Rs_ms keep m=0 for eta, but it will not used for G2
            #Rsms = [ self.Rc/(Np**(m/Np)) for m in np.arange(Np+1)]     # Rc locates upside
            Rsms = [ Np**(-m/Np) for m in np.arange(Np+1)]               # without Rc
            Rsms_1 = Rsms[1:]
            ### m+1 -> m: eta_sm = (R_sn-m - R_sn-m-1 )^2
            etasms=[ (Rsms[m+1]-Rsms[m])**(-2) for m in np.arange(Np)]
            Rs.extend(Rsms_1)
            etas.extend(etasms)
            #print(etas, Rs)
        else:
            print(f"{self.func_param} is not available")
            sys.exit(100)
        for element in elements:
            _G = make_symmetry_functions(sftype='G2', etas=etas, Rsms=Rs, elements=elements)       # Rc is out of eta
            #_G += make_symmetry_functions(type='G2b', etas=etasms, elements=elements)   # Rc is in eta
            _G += make_symmetry_functions(sftype='G4', etas=[0.005], zetas=[1., 4.], gammas=[+1., -1.], elements=elements)
            G[element] = _G
        return G

class ZN_param:
    name='zn'

class BS_param:
    name='bs'


