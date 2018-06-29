import os
from atom_valence import *       # fatom_valence
#from mplt_mo_ini import *
from common import *
from mplt_mo_ini import *
import numpy as np
import operator
import heapq
import collections
#print "\n".join(sys.path)


"""
if atom input exists, it becomes selective in MO drawing
    if natom != 1, mo_type should be input
if not atom input, default is S='all'
    if S != 'all', atom list, mo_type should be coded
"""
# Atom id such as C-1, C-2 are not coded yet
#FL_Atom_types=[             
#    [['Ni','d']],           # 1st type
#    [['C','p'],['O','p']]   # 2nd type
#    ]

def atom_compress(atoms):
    new_atoms=[]
    for atom in atoms:
        if not atom in new_atoms:
            new_atoms.append(atom)
    return new_atoms

def atom_decompose(latoms):
    latom_name=[]
    latom_index=[]
    for atom in latoms:
        m = re.search('\d', atom)       # returns matching object
        if m:
            atom_name=atom[:m.start()]
            atom_index=atom[m.start():]
            latom_name.append(atom_name)
            latom_index.append(atom_index)
        else:
            latom_name.append(atom)
            latom_index.append("")
    print latom_name, latom_index, "in", whereami()
    atom_species = atom_compress(latom_name)
    return latom_name, latom_index, atom_species

def All_atom_count(latoms):
    all_atom_count={}
    for atom in latoms:
        m = re.search('\d', atom)       # returns matching object
        if m:
            atom_name=atom[:m.start()]
            atom_index=int(atom[m.start():])
            if atom_name in all_atom_count:
                if all_atom_count[atom_name] < atom_index:
                    all_atom_count[atom_name] = atom_index
            else:
                all_atom_count[atom_name] = atom_index
        else:
            all_atom_count[atom] = 1

    print all_atom_count, "in", whereami()
    return all_atom_count

def Cal_Ncore(dict):
    ncore = 0
    for atom in dict:
        nc = Atom_Core_631g[Atom_Table.index(atom)]
        ncore += nc * dict[atom]
    print "all the number of core basis:", ncore, "in", whereami()
    return ncore        

### MO ENE block analysis
def f_extract_ene(index_homo_line, line_ene_ab):
    """
        make list of MO_ene, MO_id as for one file
        index_homo_line: line index for homo - last line of occupied enegy
        line_ene_ab: [ [Alpha: line0, line1, ... ], [Beta: line0, line1, ...]}
    """
    #### 2d means alpha, beta
    if len(index_homo_line) == 1:
        beta = 0
    else:
        beta = 1
    #print len(index_homo_line), len(line_ene_ab)
    f_imo_line2d=[]
    iab=0 
    ### 2D variables are here
    list_ene_ab=[]   # 2d for [[alpha], [beta]]
    e_homo_ab=[]
    e_homo_id_ab=[]
    ### split each lines of alpha beta
    for ab_line_list in line_ene_ab:     # alpha beta [ A[line...],B[line...]]
        # 1d variables are here
        list_ene=[]
        j=0     # ene-line index
        ihomo=0
        ### split alpha into each energy line
        for e_line in ab_line_list:     # line index, it might be list[0] / [ list[0], list[1]]
            ene=re.split("\s+", e_line)
            one_line_ene=[]             # for calculation of not full energy list in a line
            # split energy line into each energy value
            for x in ene:
                if x :                   # remove empty entry
                    list_ene.append(x)  # store in 1D
                    one_line_ene.append(x)
            ### for H-atom
            if index_homo_line[iab] == -1:
                ehomo = -10
                ihomo = 0
                e_homo_ab.append(ehomo)
                e_homo_id_ab.append(ihomo)
                continue
            #### finding homo ene level
            elif j == index_homo_line[iab]:
                ehomo=list_ene[len(list_ene)-1]  # the last energy of last occ line
                ihomo=j*Nene_line_QC+len(one_line_ene)
                e_homo_ab.append(ehomo)
                e_homo_id_ab.append(ihomo)
            j+=1
        # after finish alpha/beta            
        list_ene_ab.append(list_ene)     # store ie 2D
        iab += 1
    # after finish a/b list        
    print "HOMO energy, ID, beta, max_lumo in:", whereami(), e_homo_ab, e_homo_id_ab, beta
    return list_ene_ab, e_homo_ab, e_homo_id_ab


### MO Coefficient Block

def make_atom_basis_id(lsymbol):
    id_atom=lsymbol[0]
    atom_name=lsymbol[0]
    if re.match("\d",lsymbol[1]):
        id_atom+=lsymbol[1]        # id_atom == Ni1, Ni2, H30, H31, etc
        l_basis=lsymbol[2:]
    else:
        l_basis=lsymbol[1:]
    id_basis=l_basis[0]
    if len(l_basis) == 2:
        id_basis+=l_basis[1]
    return atom_name, id_atom, id_basis


## not complete
def moc_line_analysis(line, nenergy, crit,imo, imo_ene, p_atom_name, icore,ncore):

    moc_block=[]

    lmoc_line=line.split()
    if V_print_moc >=2: print lmoc_line
    ncol = len(lmoc_line)
    ind = ncol - nenergy
    lcoeff = lmoc_line[ind:]
    l_left = lmoc_line[0:ind]
    ibasis = l_left.pop(0)
    
    ### make atom id
    if V_print_moc >= 1: print ibasis, l_left, lcoeff
    #### FILTER 1  skip core orbital based on atom species in lmoc_line[1]: follow QChem basis
    #   as for new atom, skip number of core basis and count valence basis
    if atom_name != p_atom_name:
        p_atom_name = atom_name
        icore=0 #print lmoc_line
        ncore=Atom_Core_631g[Atom_Table.index(lmoc_line[1])]
        #print ncore
        #### if there is beta electron, makes error
        if V_print_moc >= 2: print "atom ID & num of core atomic orbital: ", id_atom, ncore
    if icore < ncore:
        icore+=1
        ### ruturn with increase imoc, skip in main
        return moc_block, p_atom_name, icore, ncore
    else:
        #### for each MOC line, SAVE nenergy==6 (if last, nenergy<=6) in moc_block
        for j in range(nenergy):
            mo_symbol=[]
            moc_tmp=[]
            #### FILTER 2 of MOC line by MOcoeff_crit
            if abs(float(lcoeff[j])) > crit:
                #mo_symbol.extend((l_left[0], l_left[1]))
                mo_symbol.extend((id_atom, id_basis))
                moc_tmp.extend((mo_symbol,imo[j],imo_ene[j],lcoeff[j]))
                                
                #### Add moc line in moc_block of 6 mo_id
                moc_block.append(moc_tmp)
                print(np.matrix(moc_block))
        ### return with the same imoc, use moc_tmp
        return moc_block, p_atom_name, icore, ncore





#   MO types and MO atoms
#   not used
def get_mo_types(i_type):   
    #   for everythin
    if i_type == 0:
        MO_atoms=[]
        mo_type="None"
    #   for CO2
    elif i_type == 1:
        MO_atoms=['C', 'O1', 'O2']
        mo_type="All"                      # print all the atoms
            
    #   for 1-PP, 2-PPP,
    elif i_type == 2:
        MO_atoms=['Ni','P']                 # draw only for Ni
        mo_type='1sub'

    #   for 3-Me-PNP, 4-PNP
    elif i_type == 3:
        MO_atoms=['Ni','P','N']
        mo_type='1sub'

    #   for CO2 complex
    elif i_type == 4:
        MO_atoms=['Ni', 'C', 'O']           # 1sub-a for additional atom id of C-1
        mo_type="Two"                    # Ni is necessary and 1 of C or O

    else:
        print "add more MO type list manually"
        exit(20)
    return  MO_atoms, mo_type    

### MO types and MO atoms
def get_atom2motype(atoms):

    if not atoms:
        mo_type = 'ALL'
    #elif len(atoms) == 1:
    #    mo_type = 'One'
    else:
        mo_type = 'SEL'

    if len(atoms) == 2:
        if set(atoms) == set(['Ni','P']):
            mo_type = 'SUB'
    elif len(atoms) == 3:
        if set(atoms) == set(['C', 'O', 'O']):
            mo_type="ALL"                      # print all the atoms
        elif set(atoms) ==set(['Ni','P','N']):
            mo_type='SUB'
    elif len(atoms) == 4:
        if set(atoms) == set(['Ni', 'C', 'O','O']):           # 1sub-a for additional atom id of C-1
            mo_type="SEL"                                   # Ni or C or O or O
    return mo_type    

def get_mo_labels(files, atoms, imotypes):
    Fl_MO_atoms=[]
    Fl_MO_type=[]
    for index in imotypes:
        #### if atoms are input, draw method is selective
        if atoms:
            if len(atoms) == 1:
                mo_type='ONE'
            else:
                #### if there are atom_list, and more than one atom, mo_type should be input
                print "Usage:: A='atom1 atom2' T[ype]= all[1sub]"
                exit(11)
        mo_atoms, mo_type = get_mo_types(index) 
        Fl_MO_atoms.append(mo_atoms)
        Fl_MO_type.append(mo_type)

    if atoms:
        print  "mo_type list: ", MO_type_list
        for type in Fl_MO_type:
            if not type in MO_type_list:
                print "type should be one of ", MO_type_list
                exit(10)
        print "atoms for MO in ", os.path.basename(__file__), Fl_MO_atoms
        print "MO_type in ", os.path.basename(__file__) ,": ", Fl_MO_type

    print "input files in module ", os.path.basename(__file__),  ": ", files

    '''
    1                                       88 homo    
    2   92(w) 93(w) 94(s,bo,O1) 99(s,bb,C-O1)           100(s,bb-ab,C-O1
    3   CO2                                             lumo
    Link between MOs are given manually
    indices: homo-1=-1 homo=0, lumo=1, lumo-1=2
    MO_link_hl_id=[[id_1st_mol, id_2nd_mol],[id_2nd_mol, id_3rd_mol]]
        1st link between 1-2 mols 2nd link between 2-3 mols
    '''
    return Fl_MO_atoms, Fl_MO_type

def get_homo_ind_in_qcdic(qc_dic, ihomo):
    ind = -1
    print qc_dic
    for i, qc_key in enumerate(qc_dic):         # default .keys()
        print i, qc_key
        if qc_key == ihomo:
            ind = i
            break
    if ind == -1:
        print "there is no homo id in tuple", whereami()
        exit(34)
    else:
        return ind
def get_homo_ind_in_qclist(qc_lclass, ihomo):
    ind = -1
    for i, qc_class in enumerate(qc_lclass):         # default .keys()
        if V_print_link >= 3: print i, qc_class.imo, whereami()
        if qc_class.imo == ihomo:
            ind = i
            break
    if ind == -1:
        print "there is no homo id in tuple", whereami()
        exit(34)
    else:
        return ind

def matched_key_bmol(sorted_bases, key1):
    """ 
        sorted bases [('dxy': 0.99), ...
        return order, coeff
    """
    coeff = 0.0
    #print sorted_bases
    for i, bas_co in enumerate(sorted_bases):
        print i, bas_co, key1, whereami()

        if bas_co[0] == key1:
            coeff = bas_co[1]
            return i+1, coeff
    return i+1, coeff

def find_level_nbases(n, qc_list, ihomo,imo_homo, tag_hl, matom):
    level_base=[]
    if tag_hl == 'homo':
        locc = qc_list[:ihomo+1]
        locc.reverse()
    elif tag_hl == 'lumo':
        locc = qc_list[ihomo:]

    for qc in locc:     # homo homo-1 homo-2 ...
        buffer = heapq.nlargest(n, qc.bas_dic.items(), lambda i: i[1])
        for mbasis in buffer:           # 1d
            if re.search(matom, mbasis[0]):
                ind = qc.imo - imo_homo                 
                level_base.append([ind,mbasis[0]])
                if tag_hl == 'homo':
                    j = ind
                elif tag_hl == 'lumo':
                    j = ind + 1
                print "%s%2d level [1st base:coeff] = %-7s %5.2f" % (tag_hl, j, mbasis[0], mbasis[1]) 
                break
    return level_base

###  MO_link_hl_id: pair between 1-2 and 2-3 fragments
def get_link(qc_lclass_A, ihomo_A, qc_lclass_B, ihomo_B, ltype, nlinking, atom_groups):
    """ 
        qc_dic = { imo: class QC_imo}
        QC_imo (self.energy, self.bc_dic
        get highest coeff basis
    """
    tag=0
    mo_hl=[]
    ### get homo index in the list: this is ordered dictionary
    h_index_A = get_homo_ind_in_qclist(qc_lclass_A, ihomo_A)
    h_index_B = get_homo_ind_in_qclist(qc_lclass_B, ihomo_B)
    print atom_groups
    ### for single atom Ni, first mol is reactant, left-side
    print "Left Mol OCC"
    AObases=find_level_nbases(1, qc_lclass_A, h_index_A, ihomo_A, 'homo', atom_groups[0])
    print "Left Mol VAL"
    AVbases=find_level_nbases(1, qc_lclass_A, h_index_A, ihomo_A, 'lumo', atom_groups[0])
    
    if ltype == 'flow':
        BObases=find_level_nbases(1, qc_lclass_B, h_index_B,ihomo_B, 'homo', atom_groups[0]) 
        BVbases=find_level_nbases(1, qc_lclass_B, h_index_B,ihomo_B, 'lumo',atom_groups[0]) 
    ### for complex homo , second mol is complex
    elif ltype == 'lsplit':
        print "Right Mol OCC"
        BObases=find_level_nbases(5, qc_lclass_B, h_index_B, ihomo_B,'homo', atom_groups[0])
        print "Right Mol VAL"
        BVbases=find_level_nbases(1, qc_lclass_B, h_index_B, ihomo_B,'lumo', atom_groups[0])
        
    ### Find matched indices for A and B molecule using index i, j [0, 1, 2.. for occ
    occ_match=[]
    val_match=[]
    for i, AObase in AObases:
        for j, BObase in BObases:
            print AObase, BObase
            if AObase == BObase:
                if occ_match and j in [ y for x, y in occ_match ]:
                    print "already %d pair exists in [%d, %d] " % (j, x, y)
                    continue
                print " ---  matched in B:",j,"-th occ MO",AObase
                pair = [i, j]
                occ_match.append(pair)
                break
        if BVbases:
            for j, BVbase in BVbases:
                if AObase == BVbase:
                    if val_match and j in [ y for x, y in val_match ]:
                        print "already %d pair exists in [%d, %d] " % (j, x, y)
                        continue
                    print " ---  matched in B lumo:",j,"-th occ MO",AObase
                    pair = [i, j]
                    val_match.append(pair)
                    break
    print occ_match, val_match
    occ_match.extend(val_match)
    return occ_match
    """
    if ltype == 'lsplit':
        ### scan occ of left molecule: from homo downward

        while 1:
            if qc_class_b[j].basis == basis1:
                match_j = j
                break
            j -= 1
        homo_pair=[qc_class_a[ihomo_a].imo, qc_class_b[match_j].imo]
        mp_homo_a = 0
        mp_homo_ab=qc_class_b[match_j].imo - homo_id_b
        mo_hl.append([0, mp_homo_ab])
        ### scan lumo upward in list 2
        k= ihomo_b+1
        while 1:
            if qc_class_b[k].basis == basis1L:
                match_k = k
                break
            k +=1
        lumo_pair=[qc_class_a[ihomo_a+1].imo, qc_class_b[match_k].imo]
        mp_lumo_a = 1
        mp_lumo_ab = qc_class_b[match_k].imo - homo_id_b
        mo_hl.append([1, mp_lumo_ab])
        pass
    elif ltype == 'flow':
        ### scan homo downward in list 2
        j =  ihomo_b
        while 1:
            if qc_class_b[j].basis == basis1:
                match_j = j
                break
            j -= 1
        homo_pair=[qc_class_a[ihomo_a].imo, qc_class_b[match_j].imo]
        mp_homo_a = 0
        mp_homo_ab=qc_class_b[match_j].imo - homo_id_b
        mo_hl.append([0, mp_homo_ab])
        ### scan lumo upward in list 2
        k= ihomo_b+1
        while 1:
            if qc_class_b[k].basis == basis1L:
                match_k = k
                break
            k +=1
        lumo_pair=[qc_class_a[ihomo_a+1].imo, qc_class_b[match_k].imo]
        mp_lumo_a = 1
        mp_lumo_ab = qc_class_b[match_k].imo - homo_id_b
        mo_hl.append([1, mp_lumo_ab])
        print mo_hl

    """


###  MO_link_hl_id: pair between 1-2 and 2-3 fragments
def get_link_QCIMO(qc_class_a, homo_id_a, qc_class_b, homo_id_b, ltype, nlinking, atom_groups):
    """ use imo index """
    tag=0
    mo_hl=[]
    # if not sorted
    #qc_class_a.sort(key=operator.attrgetter('imo')
    #qc_class_a.sort(key=operator.attrgetter('imo')
    # if sorted
    ### find index of homo in the class list 1 & 2
    for i in range(len(qc_class_a)):
        if qc_class_a[i].imo == homo_id_a:
            ihomo_a = i
            print ihomo_a
            break
    basis1 = qc_class_a[i].basis
    basis1L = qc_class_a[ihomo_a+1].basis
    for i in range(len(qc_class_b)):
        if qc_class_b[i].imo == homo_id_b:
            ihomo_b = i
            print ihomo_b
            break

    if ltype == 'flow':
        ### scan homo downward in list 2
        j =  ihomo_b
        while 1:
            if qc_class_b[j].basis == basis1:
                match_j = j
                break
            j -= 1
        homo_pair=[qc_class_a[ihomo_a].imo, qc_class_b[match_j].imo]
        mp_homo_a = 0
        mp_homo_ab=qc_class_b[match_j].imo - homo_id_b
        mo_hl.append([0, mp_homo_ab])
        ### scan lumo upward in list 2
        k= ihomo_b+1
        while 1:
            if qc_class_b[k].basis == basis1L:
                match_k = k
                break
            k +=1
        lumo_pair=[qc_class_a[ihomo_a+1].imo, qc_class_b[match_k].imo]
        mp_lumo_a = 1
        mp_lumo_ab = qc_class_b[match_k].imo - homo_id_b
        mo_hl.append([1, mp_lumo_ab])
        print mo_hl

    elif ltype == 'split':
        ### scan 
        j =  ihomo_b
        while 1:
            if qc_class_b[j].basis == basis1:
                match_j = j
                break
            j -= 1
        homo_pair=[qc_class_a[ihomo_a].imo, qc_class_b[match_j].imo]
        mp_homo_a = 0
        mp_homo_ab=qc_class_b[match_j].imo - homo_id_b
        mo_hl.append([0, mp_homo_ab])
        ### scan lumo upward in list 2
        k= ihomo_b+1
        while 1:
            if qc_class_b[k].basis == basis1L:
                match_k = k
                break
            k +=1
        lumo_pair=[qc_class_a[ihomo_a+1].imo, qc_class_b[match_k].imo]
        mp_lumo_a = 1
        mp_lumo_ab = qc_class_b[match_k].imo - homo_id_b
        mo_hl.append([1, mp_lumo_ab])
        pass


    else:
        print "link type error in", whereami

    #test = [[0,0], [1,1]]
    return mo_hl

### not used
#   Fl_nlumo: file_list for maximum num of lumo for each file
def get_nlumo_linkid(nfiles):
    
    if nfiles == 1:
        #x_max=1
        Fl_nlumo=[Nmax_virtual]
        MO_link_hl_id=[]
    #### make link if nfiles >= 2    
    elif nfiles == 2:
        #x_max=2
        Fl_nlumo=[Nmax_virtual,Nmax_virtual]
        MO_link_hl_id=[[[0,0],[1,1]]]    # does it need for 2 files: CO2 vs bent CO2
    elif nfiles == 3:
        #x_max=3
        if re.search("1-P", files[0]): 
            #MO_link_hl_id=[[[0,1],[0,0],[-1,-1],[-2,-2],[-3,-3],[-4,-4],[-2,-6]], [[1,1],[0,1],[-1,0],[-2,-1],[-3,-1],[-4,0],[-6,-1],[-32,0]]]    # 1-PP w. CO2-bended
            MO_link_hl_id=[[[0,2],[0,0],[-1,-1],[-2,-2],[-3,-3],[-4,-4]], [[2,1],[0,1],[-1,0],[-2,-1],[-3,-1],[-4,0]]]    # 1-PP G631s
            #MO_link_hl_id=[[[0,2],[0,0]], [[2,1],[0,1]]]    # 1-PP G631s
            Fl_nlumo=[1,2,2]        # cut number of lumo but include degeneracy of lumo of CO2
        elif re.search("2-P", files[0]): 
            MO_link_hl_id=[[[0,0],[0,5]],[[0,1],[5,1]]]
            Fl_nlumo=[2,5,2]        # cut number of lumo but include degeneracy of lumo of CO2
        elif re.search("3-P", files[0]):
            MO_link_hl_id=[[[0,0],[0,6]],[[0,1],[6,1]]]
            Fl_nlumo=[2,6,2]        # cut number of lumo but include degeneracy of lumo of CO2
        else:
            Fl_nlumo=[Nmax_virtual,Nmax_virtual,Nmax_virtual]
    elif nfiles == 4:
        #x_max=4
        print "error"
        exit(10)
    elif nfiles == 5:
        #x_max=5
        if re.search("1-P", files[0]): 
            MO_link_hl_id=[[[0,0],[1,1]],[[0,2],[0,0],[-1,-1],[-2,-2],[-3,-3],[-4,-4]], [[2,1],[0,1],[-1,0],[-2,-1],[-3,-1],[-4,0]],[[-1,0],[0,0],[1,1],[2,1]]]    # 1-PP G631s
            Fl_nlumo=[1,1,2,2,2]        # cut number of lumo but include degeneracy of lumo of CO2
        elif re.search("2-P", files[0]): 
            MO_link_hl_id=[[[0,0],[1,1]],[[0,0],[0,5]],[[0,1],[5,1]],[[-1,0],[0,0],[1,1],[2,1]]]
            Fl_nlumo=[2,2,5,2,2]        # cut number of lumo but include degeneracy of lumo of CO2
        elif re.search("3-P", files[0]):
            MO_link_hl_id=[[[0,0],[1,1]],[[0,0],[0,6]],[[0,1],[6,1]],[[-1,0],[0,0],[1,1],[2,1]]]
            Fl_nlumo=[2,2,6,2,2]        # cut number of lumo but include degeneracy of lumo of CO2
        else:
            Fl_nlumo=[Nmax_virtual,Nmax_virtual,Nmax_virtual,Nmax_virtual,Nmax_virtual]
    else:
        print "Error: too many files of ", nfiles, " : ", os.path.basename(__file__)
        exit(11)
    #print FL_cut_nlumo
    #exit(99)
    return Fl_nlumo, MO_link_hl_id



