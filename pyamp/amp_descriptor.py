### train descriptor
### Gaussian Descriptor

import numpy as np
from amp.descriptor.gaussian import make_symmetry_functions
import sys

class GS_param:
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
        #etas = np.logspace(np.log(0.5), np.log(5), num=6)
        if self.func_param == 'log10':
            etas = np.logspace(np.log10(self.pmin), np.log10(self.pmax), num=self.nparam)
        else:
            print(f"only {scale} is available")
            sys.exit(100)
        for element in elements:
            _G = make_symmetry_functions(type='G2', etas=etas, elements=elements)
            _G += make_symmetry_functions(type='G4', etas=[0.005], zetas=[1., 4.], gammas=[+1., -1.], elements=elements)
            G[element] = _G
        return G

class ZN_param:
    name='zn'

class BS_param:
    name='bs'


def make_Gs(image, f_param='log10',  min=0.05, max=5, num=4):
    """
        Amp defaults: scale=log10, min=0.05, max=5, num=4
    """
    ### obtain atom symbols
    elements = set([atom.symbol for atom in image])
    elements = sorted(elements)     # because set is not subscriptable
    G = {}
    #etas = np.logspace(np.log(0.5), np.log(5), num=6)
    if f_param == 'log10':
        etas = np.logspace(np.log10(min), np.log10(max), num=num)
    else:
        print(f"only {scale} is available")
        sys.exit(100)
    for element in elements:
        _G = make_symmetry_functions(type='G2', etas=etas, elements=elements)
        _G += make_symmetry_functions(type='G4', etas=[0.005], zetas=[1., 4.], gammas=[+1., -1.], elements=elements)
        G[element] = _G
    return G



