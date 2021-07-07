### train descriptor
### Gaussian Descriptor

import numpy as np
import sys
from common import whereami
#from amp_gsversion import g_version
from amp.descriptor.amp_source_check import Lflow, nflow
#if g_version == 2:
#    from amp.descriptor.gaussian2 import make_symmetry_functions
#elif g_version == 3:
#    from amp.descriptor.gaussian3 import make_symmetry_functions

class GS_param:
    '''
    func_param: 'log10', powNN=N^(m/N)
    Now eta & Rsm is calculated without Rc
    '''
    name='gs'
    def __init__(self, func_param='log10', pmin=0.05, pmax=5, pnum=4, cutoff = 6.5):
        self.func_param = func_param
        self.pmin = pmin
        self.pmax = pmax
        self.nparam = pnum
        self.cutoff = cutoff

    def make_Gs(self, image):
        ### obtain atom symbols
        elements = set([atom.symbol for atom in image])
        elements = sorted(elements)     # because set is not subscriptable
        G = {}
        Np = self.nparam
        Rc = self.cutoff
        #etas = np.logspace(np.log(0.5), np.log(5), num=6)
        if self.func_param == 'log10':
            print(f"parameter function {self.func_param} was chosen")
            etas = np.logspace(np.log10(self.pmin), np.log10(self.pmax), num=Np)
        elif self.func_param == "powNN":
            ### G2a: for N, num(m) = N+1, m={0,1,...,N)
            ### exclude Rc and use gaussian2.py
            #if g_version == 2:
            print(f"g_version {g_version}: Rc stays in gaussian2")
            ### etas with Rs = 0
            etas = [ Np**(2*m/Np) for m in np.arange(Np+1)]            # without Rc
            Rs = list(np.zeros(len(etas)))
            ### etas with Rs != 0; Rs_ms keep m=0 for eta, but it will not used for G2
            Rsms = [ Np**(-m/Np) for m in np.arange(Np+1)]               # without Rc
            Rsms_1 = Rsms[1:]
            etasms=[ (Rsms[m+1]-Rsms[m])**(-2) for m in np.arange(Np)]
            ### m+1 -> m: eta_sm = (R_sn-m - R_sn-m-1 )^2
            ### Bug exists: Rc is included here for gaussian3.py
            #elif g_version == 3:
             #   print(f"g_version {g_version}: Rc is included in etas and Rs")
             #   etas = [ (Np**(m/Np)/Rc)**2 for m in np.arange(Np+1)]           
              #  Rs = list(np.zeros(len(etas)))
              #  Rsms = [ Rc/(Np**(m/Np)) for m in np.arange(Np+1)]       # Rc locates upside
               # Rsms_1 = Rsms[1:]
               # etasms=[ (Rsms[m+1]-Rsms[m])**(-2) for m in np.arange(Np)]
            #else:
            #    print(f"module import error in {__name__}")
            #    sys.exit(0)
            Rs.extend(Rsms_1)
            etas.extend(etasms)
            print(etas, Rs, f"in {whereami()} in {__name__}")
        else:
            print(f"{self.func_param} is not available")
            sys.exit(100)
        for element in elements:
            _G = make_symmetry_functions(sftype='G2', etas=etas, Rsms=Rs, elements=elements)       # Rc is out of eta
            _G += make_symmetry_functions(sftype='G4', etas=[0.005], zetas=[1., 4.], gammas=[+1., -1.], elements=elements)
            G[element] = _G
        nflow.ncount += 1
        if Lflow: print(f"finished making Gs {nflow.ncount}: {whereami()}() in module {__name__}")
        return G

class ZN_param:
    name='zn'

class BS_param:
    name='bs'


