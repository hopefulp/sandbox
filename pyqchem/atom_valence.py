import re

Atom_Table      =[ 'H', 'B', 'C', 'N', 'O', 'F', 'P', 'Ni' ,'Fe']
#Atom_Valence=[ 2,  2,   2,   2,   3,   4 ]
Atom_Core_631g  =[  0,   1,   1,   1,   1,   1,   5,   9, 9 ]    
######## 6-31Gs
######## Ni: 9
#### P : 1s 2s 2p*3         = 5
#### Ni: 1s 2s 2p*3 3s 3p*3 = 9


Atom_checked=[]
Atom_val_dic={}
Atom_Core_dic={}
#### obtain atoms from filename
def fatom_valence(fname):
    for atom in Atom_Table:
        tag_atom=0
        if re.search(atom, fname):
            #### remove atom identity of B for -B- in filename
            nonatom='-'+atom+'-'
            if re.search(nonatom, fname):
                tag_atom=1
            if tag_atom == 0:                    
                #print atom
                Atom_checked.append(atom)
                Atom_Core_dic[atom]=Atom_Core_631g[Atom_Table.index(atom)]
    #return Atom_checked
    return Atom_Core_dic


