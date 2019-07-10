import numpy as np

"""
    quantum calculation
    class MO-level: basis - coefficient
"""

class AtomBasis_index:
    """ given basis index: have atom_name, id_atom, id_basis """
    def __init__(self, ibasis, atom_name, atom_name_index, basis):
        self.index_basis=ibasis
        self.atom_name=atom_name
        self.atom_name_id=atom_name_index
        self.basis=basis
        self.mo_symbol=(self.atom_name_id, self.basis)

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
        self.basis = []     # basis, coeff is not used now
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

def dump_moc(l_moc):
    list_2d=np.array(l_moc)
    print("Rank", list_2d.shape)
    for moc_block_ele in l_moc:
        print (moc_block_ele)


