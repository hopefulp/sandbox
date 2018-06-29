"""
    quantum calculation
    class MO-level: basis - coefficient
"""

class QC_IMO_v1:
    """ contains qmo id and atomic basis - coeff as dictionary """
    def __init__(self, nid, basis, coeff, energy):
        self.imo = nid
        self.basis = basis
        self.coeff = coeff
        self.ene = energy

class QC_IMO:
    """ contains qmo id and atomic basis - coeff as dictionary """
    ### initialize by imo & energy
    def __init__(self, nid, energy, ndic):
        self.imo = nid
        self.energy = energy
        ### list in the order: if dictionary, treat the same keys
        self.basis = []
        self.coeff = []
        self.bas_dic = ndic
        
    def adder_list(basis, coeff):
        self.basis.append(basis)
        self.coeff.append(coeff)

    def adder_dict(basis, coeff):
        self.bas_dic[basis] = coeff

    def input_dict(dict):
        self.bas_dic=dict

### save imo to dictionary w. keys of imo & values of class of basis-coeff pairs
    
class QC_imo:
    """ contains qmo id and atomic basis - coeff as dictionary """
    ### initialize by imo & energy
    def __init__(self, energy, dic):
        self.energy = energy
        ### list in the order: if dictionary, treat the same keys
        self.basis = []
        self.coeff = []
        self.bc_dic = dic
        
    def adder_list(basis, coeff):
        self.basis.append(basis)
        self.coeff.append(coeff)

    def adder_dict(basis, coeff):
        self.bas_dic[basis] = coeff

    def input_dict(dict):
        self.bas_dic=dict



