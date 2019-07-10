#### Module to move block of main script to here
import re
import os
import sys
from atom_valence import *
from common import *
from mplt_mo_ini import *
from mp_mo_dump import *
from qmo import *

### keywords for MO Energies
# first "Occupied" for Alpha, second for beta
KW_MOene_occ="Occupied"
KW_MOene_vir="Virtual"
KW_MOene_beta="Beta MOs"
KW_MOene_end="-----"
#### key word
#key_MOa="Alpha MOs"
#key_MOb="Beta MOs"

################################# BLOCK 1 :: MO energy##############################################
### obtain MO energy lines, then
### extract energies, homo and their imo index

def get_MO_energies(f):
    """
    input file pointer
    return energies, indices, homo
    """
    flag_block = "OFF"
    eline_alpha=[]
    eline_beta=[]
    i_hline_ab=[]

    for line in f.readlines():
        ### Find MO block & initialize (first "occupied" is for Alpha spin)
        if flag_block == "OFF":
            if re.search(KW_MOene_occ, line):
                print ("####################  BLOCK 1 (MOE):: obtain MO Energy in {}:{}()".format(os.path.basename(__file__), whereami()))
                flag_block="ON"
                raise_beta=0              # beta option is initiated
                ### initiate all the variables used in this block
                i_eline=0            # energy line index starts from 0
                tag_beta="OFF"
                flag_vline="NO"         # now occupied orbital
                continue
        #### inside MO Energy Block
        elif flag_block == "ON":
            ### if options except energies, find keyword in reverse order (?)
            ### 1. find end block line: if arrived the END line of MO energy
            if re.search(KW_MOene_end, line):
                break               # end of Block
            #### 2. skip non-ene line (symmetry lines: watch out B|Beta)
            elif re.search("A", line): # symmetry simbols: A1, B1, A', A''
                continue
            ###### 2nd Alpha Virtual Block 
            #### 3. find Virtual line
            elif re.search(KW_MOene_vir, line):
                flag_vline="YES"                            # Virtual orbital energies
                #   for H atom (one electron - no beta occ)
                #if raise_tag == 1:
                #    if not f_line_ene_b:
                #        print (f_line_ene_b, "is false")
                #        i_hline=-1
                ### at Virtual line: add homo line index
                i_hline_ab.append(i_hline)                  # twiced if beta_tag, i_hline is assigend below
                if vprint_moe >=1: 
                    print ("HOMO line index: i_hline = ", i_hline, end=' ')
                    if tag_beta == "YES": 
                        print("beta")
                    else: 
                        print("alpha")
                continue
            #### in case of Beta, there is blank line, how to treat BLANK
            #### 4. find blank for beta
            elif not line.strip():
                flag_vline="NO"     # virtual energy ended and prepare for beta
                continue
            #### 5. find BETA
            elif re.search(KW_MOene_beta, line):
                raise_beta=1                            # turn on beta tag
                tag_beta="YES"
                i_eline=0                                   # reset energy line index, it starts from 0
                continue
            ### 6. find "Occupied" again for beta
            elif re.search(KW_MOene_occ, line):
                continue
            #### save energy line if skipped all the previous if condition
            else:
                if tag_beta == "YES":
                    ### if beta is detected but is the same as alpha, igmore for bent CO2 in 1-PP-B
                    #if ib_eline == 0:
                    #    print ("check beta line", line.rstrip, f_line_ene_a[ib_eline])
                    #    if line.rstrip() == f_line_ene_a[ib_eline]:
                    #        print ("Beta spin is detected but same as Alpha spin")
                    #        Fl_beta_tag[fid]=0
                    #        tag_MO="Done"
                    #        continue
                    eline_beta.append(line.rstrip())
                    if flag_vline != "YES":
                        i_hline=i_eline
                    #ib_eline+=1
                ### 1st block for Alpha, occupied
                else:
                    eline_alpha.append(line.rstrip())
                    if flag_vline != "YES":
                        i_hline = i_eline
                i_eline += 1
                continue
    ### Wrap up MO Block
    """
        eline_alpha=[l1, l2...]
        eline_2d =[[Alpha:l1, l2,...],[Beta:l1, l2,...]]
        FL_homo_ene =[[f1:[A-homo, B-homo],[f2:...
        Fl_homo_ene =[f1-homo, f2-homo,...]
        FL_homo_id  =...
    """
    eline_2d=[eline_alpha, eline_beta]
    #print(eline_2d, i_hline_ab)
    fL_ene_ab,f_homo_ab,f_homo_id_ab=extract_ene_4line(i_hline_ab, eline_2d)
    #print(len(fL_ene_ab), len(i_hline_ab),fL_ene_ab, i_hline_ab)
    '''
    FL_homo_ene.append(f_homo_ab)
    print (FL_homo_ene)
    if Fl_beta_tag[fid] == 1:
        Fl_homo_ene.append( max(FL_homo_ene[fid][0], FL_homo_ene[fid][1]) )
    else:
        Fl_homo_ene.append( FL_homo_ene[fid][0])
    FL_homo_id.append(f_homo_id_ab)
    FL_MOene_list.append(fL_ene_ab)
    
    f_norbital_ab.append(len(fL_ene_ab[0]))
    if Fl_beta_tag[fid] == 1:
        f_norbital_ab.append(len(fL_ene_ab[1]))
    FL_norbital.append(f_norbital_ab)
    '''
    #return eline_2d,i_hline_ab,raise_beta
    return fL_ene_ab,f_homo_ab,f_homo_id_ab,raise_beta

### MO ENE block analysis: from Energy lines extract ene and id
def extract_ene_4line(index_homo_line, line_ene_ab):
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
    print(f"{whereami()}(): HOMO energy, ID, beta \n\t\t   {e_homo_ab}  {e_homo_id_ab}  {beta}  ")
    return list_ene_ab, e_homo_ab, e_homo_id_ab

### keywords for MOC 
KW_MOcoeff="MOLECULAR ORBITAL COEFFICIENTS" # common for alpha & beta
KW_nbasis="basis functions"

#################################################### BLOCK 2 :: MO Coeff ########################################
### ene & index are obtained here
### find key_word_MOcoeff and initialize
def obtain_MOC(f,beta,homo_id_ab,norb_ab,atomlist_id=None,mo_type=None):
    """
    input:: f = file pointer
            homo_id_ab: homo index of alpha, beta of input file
            norb_ab: maximum number of orbitals of alpha, beta of input file
    return moc
    """
    global n_loop
    print (f"####################  BLOCK 2 (MOC) :: MO Coeff starts in {whereami()}() ")
    print(f"{whereami()}()::  atomlist: {atomlist_id}; mo_type: {mo_type}")
    flag_block = "OFF"
    flag_moc_beta = "OFF"
    tag_beta='OFF'

    all_id_atoms=[]
    base_dic={}
    moc_block=[]         # keeps all info in one block
    input_atoms=[]       # obtain from input atomlist
    input_atoms_id=[]
    imo_sel_ab=[]
    ene_sel_ab=[]
    ### for alpha, beta list [[alpha],[beta]]
    L_imo_select=[]     
    L_ene_select=[]
    L_bcdic_select=[]    
    l_QCIMO=[]           # l_QCIMO=[class list]
    L_QCIMO=[]           # l_QCIMO=[class list]
                         # class QC_IMO {self.imo, self.ene, self.basis=[], self.coeff=[], self.bas_dic(redundant, if convenient)

    k = 0               # imoc was converted to k: 0 ~ nbase+1, total nbase+2 lines
    i_block = 0          # imoc_block was converted to i_block
    n_check = 0
    Tag_fine = 0
    for line in f.readlines():
        ### Find nbasis for MO Coefficient in advance
        if re.search(KW_nbasis, line) and re.search("shells", line):
            field=re.split("\s+", line)
            field=[x for x in field if x]
            nbasis=float(field[5])
            #print "found num of basis: ", nbasis
            continue
    
        ###### FIND MOC block & initialize
        if flag_block == "OFF":
            ###### FOUND MOC PART & INITIALIZE
            if re.search(KW_MOcoeff, line) :    #common for alpha & beta
                if vprint_moc >=1: print("MOC found w. beta_tag \"{tag_beta}\": start w. initialization")
                flag_block="YES"
                #### constants from module mplt_mo_ini
                max_nlumo = Nmax_virtual
                #### obtain atoms for MO and mol types in QChem outfile
                if atomlist_id:
                    #atomlist_id=atomlist.split()   # atomlist="Ni C1", atomlist_id=["Ni", "C1"]
                    #print(atomlist_id, atomlist)
                    atom_style='name'
                    # Ni is 'name' 'Ni, C1' is 'index'
                    for atomid in atomlist_id:
                        if re.search('\d', atomid):
                            atom_style='index'
                    if atom_style=='name':
                        input_atoms=atomlist_id[:]
                        input_atoms_id=input_atoms    # in case no index, input_atoms == input_atoms_id
                    ### USE index
                    elif atom_style=='index':
                        latom_name, latom_index,atom_species= atom_decompose(atomlist_id)
                        input_atoms=latom_name[:]
                        input_atoms_id=atomlist_id
                if vprint_moc >=1: print(f"atom_style is '{atom_style}' for atomlist_id {atomlist_id}")

                if mo_type:
                    pass
                else:
                    mo_type=get_atom2motype(atoms)         

                #mo_type.append(f_mo_type)             #Fl_MO_type[:]
                #print (fl_mo_atoms_id, atoms_id, f_mo_type, mo_type, "in", whereami())
                ##atoms, Fl_MO_type = get_mo_labels(files, atoms, i_mo_types)
                ##Fl_nlumo, MO_link_hl_id = get_nlumo_linkid(nfile)
                ##print fl_mo_atoms, f_mo_type, Fl_nlumo, MO_link_hl_id
                ### use common homo index for a|b
                #fL_homo_id=FL_homo_id[fid]      # list of one|two elements
                ## Alpha has always large value: may be not
                #ilumo_cut=homo_id_ab+Fl_nlumo[fid]
                if tag_beta == 'OFF':
                    ihomo=homo_id_ab[0]
                else:
                    ihomo=homo_id_ab[1]
                ilumo_cut=ihomo + max_nlumo
                ### define lumo index cut
                if norb_ab[0] < ilumo_cut:
                    ilumo_cut = norb_ab[0]
                print (f"HOMO id with beta-tag \"{tag_beta}\" & Virtual Orbital cut id:: {ihomo} {ilumo_cut} ")
                if tag_beta == "ON":
                    f_homo_id=homo_id_ab[1]
                    print ("###################  BLOCK 2 :: MO Coeff beta starts ")
                    k=0
                    i_block=0
                    imoc_valuable=[]
                    moc_block=[]
                    imo_sel_ab=[]
                    ene_sel_ab=[]
                #print (" MO Levels")
                continue
        #################### INSIDE MOC PART:: including all the MOC blocks
        elif flag_block == "YES":
            ######### Deal With MOC PART: as unit of MOC BLOCK (nbasis+2 lines)
            ###### escape if end of MO Coefficients
            if vprint_loop>=2: print(f"k increase {k}")                       # this is OK
            if not line.strip():
                if beta == 0:
                    print ("End of MOC of alpha")
                    break
                else:
                    print("Check:: MOC beta has been read")
                    break
            ##################### READ 1ST LINE 
            ##### save MO id 
            if k == 0:                                       # k runs 0 - nbasis+1, inside MOC 1 block
                l_imo=line.split()                           # k: total nbasis + 2 lines
                Ne_col_block=len(l_imo)                      # Ne_clnm_block is num of ene column in each block, recalculate each time
                if vprint_moc >=2 : print ("k==0; imo =", line.strip())
                if vprint_log >=1 and n_check==0: print(f"line-k:{k}-{n_check}: have read index of mo, imo: {l_imo}"); n_check+=1
                #print(l_imo)
                k += 1
                ##### IN NEW Block initialize here
                p_imo=int(l_imo[0])
                imoc_valuable=[]
                continue
            ##################### READ 2ND LINE
            ##### save energy
            elif k == 1:
                ### if symmetry line, skip
                #if re.search("A1", line) or re.search("B1", line):
                if re.search("A", line):
                    continue
                imo_ene=line.split()
                imo_ene.pop(0)                              # first column is "eigenvalues:"
                if vprint_moc >= 2: print ("k==1; imo_ene", imo_ene)
                if vprint_log >=1 and n_check==1: print(f"line-k:{k}-{n_check}: have got moe energy in the block: {imo_ene}"); n_check+=1
                k += 1
                ##### IN NEW Block initialize here
                p_moe=float(imo_ene[0])
                continue
            # if's column
            ##################### READ REMAINING ALL THE LINES
            #### Analize ONE moc line for one basis in a MOC_BLOCK: from 2 to 1+nbasis line
            elif k <= 1+nbasis:                                      # from k==2, nbasis starts
                # moc_line_analysis(line, Ne_col_block, temp_atom_name)
                lmoc_line=line.split()
                if vprint_moc >=2 and k==2: print ("moc-i {}".format(k-1), lmoc_line)
                # each row: ncol is different
                ncol = len(lmoc_line)
                ind = ncol - Ne_col_block
                lcoeff = lmoc_line[ind:]
                l_symbol = lmoc_line[0:ind]
                ind_basis = l_symbol.pop(0)     # use this one to get AtomBasis_id.mo_symbol
                #### base_dic{base id: class AtomBasis_index} is made at once
                #print(f"k {k}", line)
                if i_block == 0:
                    ### make atom id: id_atom=Ni1, C1, O1, O2, H13 etc; id_basis=s, px, dxx, f1, f2 etc
                    atom_name, id_atom, id_basis = make_atom_basis_id(l_symbol)
                    #print (l_symbol, whereami(),atom_name,id_atom,id_basis)
                    tmp_base = AtomBasis_index(ind_basis, atom_name, id_atom, id_basis)
                    if not id_atom in all_id_atoms:
                        all_id_atoms.append(id_atom)
                    base_dic[ind_basis]=tmp_base
                #if vprint_moc >= 1: print ind_basis, l_symbol, lcoeff, atom_name, id_atom
                ############### MOC_BLOCK FORMAT:: mo_symbol, imo, moe, moc
                ##### each line is decomposed into 6 MOC format
                if vprint_log >=1 and n_check==2: print(f"line-k:{k}-{n_check}: decompose line to get 6 items in moc_block[]"); n_check+=1
                for j in range(Ne_col_block):
                    mo_symbol=[]
                    moc_tmp=[]
                    #### FILTER 1: 6 data of MOC line filtered by MOcoeff_crit
                    #### make moc_block for one MOC block w. 6 mo's
                    if abs(float(lcoeff[j])) > MOcoeff_crit:
                        #print(base_dic[ind_basis].mo_symbol,l_imo[j],imo_ene[j],lcoeff[j])
                        moc_tmp.extend((base_dic[ind_basis].mo_symbol,l_imo[j],imo_ene[j],lcoeff[j]))
                        ### moc_block contains required data type for 1 MOC block ( Ne_col_block*nbasis maximum)
                        if vprint_log >=1 and n_check==3:
                            print("moc format:: mo_symbol(class) mo-index, mo-energy, MOCeff")
                            print(f"\t{moc_tmp} in block scanning") 
                            print("\t'moc_block[]' gathers moc format for a block")
                        moc_block.append(moc_tmp)
                        #print(np.matrix(moc_block))
            # if inside block: k <= 1+nbasis
                ############################ HAVE READ ONE BLOCK
                ##### Have read 6 imo and saved in moc_block[], so Treat it and INITIALIZE variables again
                if k == 1+nbasis:
                    if vprint_log >=1 or vprint_loop>=1: print(f"line-k:{k}-{n_check}:{i_block}-th Block was read and is treated"); ### leave for error check 
                    ### for first MOC block, calculate core
                    if i_block == 0:
                        all_atom_count = All_atom_count(all_id_atoms)
                        ncore = Cal_Ncore(all_atom_count)
                        if vprint_moc >=1:
                            print(f"In {i_block}-th MOC block:: moc format::")
                            print(f"\t{moc_tmp} in block scanning")
                            print("\t'moc_block[]' gathers moc format for a block")

                    #### treat block of MOC here
                    #### new style:           [['Ni4', 'd4'], '100', '-0.016', '0.19560']
                    #### FORMAT of moc_block: [[atom_id,basis_name], mo_id, mo_ene, mo_coeff
                    #### SORTED by  mo_id
                    sorted_moc_block=sorted(moc_block, key=lambda x:float(x[1]))
                    if vprint_log >=1 and n_check==3: 
                        for  x in sorted_moc_block:
                            print ("in sorted moc_block: ", x)
                        n_check+=1
            # if inside block: k <= 1+nbasis
                # if k == 1+nbasis
                    if vprint_log >= 1 and n_check==4: print (f"line-k:{k}-{n_check}:treat sorted moc_block[] in for-loop") 
                    i_sortb=0
                    nmocblock=len(sorted_moc_block)
                    ################## SCAN SORTED MOC_BLOCK[] in moc-format
                    for moc_block_ele in sorted_moc_block:      # to the end
                        i_sortb+=1
                        b_imo       = int(moc_block_ele[1])           # moc_block[1] == str(imo)
                        b_atom_id   = moc_block_ele[0][0]
                        b_basis_id  = moc_block_ele[0][1]
                        b_mo_ene    = float(moc_block_ele[2])
                        b_coeff     = moc_block_ele[3]
                        b_atomid_n_basisid='-'.join(moc_block_ele[0])
                        ########## SKIP IF imo <= ncore
                        #### FILTER 2  skip lower part of core orbital based on atom species in lmoc_line[1]: follow QChem basis
                        if b_imo <= ncore:            # NCORE is treated here
                            if b_imo > p_imo:
                                if vprint_log >=1: print(f"{b_imo} imo is skipped due to core")
                                p_imo=b_imo
                                p_moe=b_mo_ene
                            continue
                        elif b_imo == ncore:
                            p_imo=b_imo
                            p_moe=b_mo_ene
                            
                        ########## FINISH IF imo > ilumo_cut
                        #### FILTER 3: cut upper part than ilumo_cut
                        #### This is for last Block w.r.t. lumo cut, when reached lumo_cut, save and finish
                        if b_imo > ilumo_cut:
                            #### dump lastly MO Coefficient before break
                            if vprint_log >=1: print(f"last adding moc block {b_imo-1} append w. beta_tag {tag_beta}") #, l_imo_sel_ab)
                            #### save imo & ene for alpha or beta
                            L_imo_select.append(imo_sel_ab)
                            L_ene_select.append(ene_sel_ab)
                            L_bcdic_select.append(ene_sel_ab)
                            L_QCIMO.append(l_QCIMO)
                            if vprint_loop>=1:print(f"rank check for L_imo_select: when b_imo > ilumo_cut, {len(L_imo_select)}")
                            Tag_fine = "OK"
                            break
                            #print("check: ", L_bcdic_select)
                            ### Finalize MOC
                            ### Beta is on construction
                            #print fL_imo_select[0]
                            #if Fl_beta_tag[fid]==0:
                            #    print ("Cut by nlumo")
                            #elif Fl_beta_tag[fid]==1 and flag_moc_beta=="ON":
                            #    break
                            #### run beta again
                            #elif Fl_beta_tag[fid]==1 and flag_moc_beta=="OFF":
                            #    flag_moc_beta="ON"
                            #    flag_moc="OFF"       # find KW for beta
                            #### FINISH MO Coeff selection when b_imo == ilumo_cut
                            ### be careful:: here is for-loop in sorted_moc_block[] and will loop continue for the remaining Block
                            #break
                        ######## DUMP for selected atom line in selected 6 mo-id's ##########
                        ####### FILTER 4       :: cut by atom species and ID here
                        ####    f_imo_dump()   :: decision for dumping the imo-block
                        ####    moc_block_ele== [['Ni4', 'd4'], '100', '-0.016', '0.19560']
                        #### 'pass' writes the line in dump string      of str_dump
                        #### 'continue' skips the line before writing   in str_dump
                    # for sorted_moc_block[]
                        #? MO type is necessary?
                        if not mo_type == None:
                            if mo_type == 'ALL':
                                pass
                            elif mo_type == "ONE":
                                if b_atom_id == input_atoms_id[0]:
                                    if vprint_filter >= 2: print (b_atom_id, atoms_id[0], "in filter 4")
                                else:
                                    pass
                            #### 'ALL' and '1sub' works the same but different in f_imo_dump()
                            #### more detail screening difference between 'ALL' and '1sub', refer to imo_dump
                            #### for more options for atom id w. index
                            elif (re.search("sel",mo_type,re.IGNORECASE) or re.search("sub",mo_type, re.IGNORECASE)):
                                if b_atom_id in input_atoms_id:
                                    pass
                                else:
                                    pass 
                            else:
                                print ("Error:: there is no atom type in", whereami() )
                                exit(11)
                        #print (moc_block_ele)
                        #### if not last ele in sorted_moc_block[]
                    # for-loop in sorted moc_block[]
                        if i_sortb != nmocblock:      # ? treat last line
                            #### index proceeds here
                            ### gather data for one imo
                            if b_imo == p_imo:
                                if float(b_coeff) < 0:
                                    dump_str="      "+b_atomid_n_basisid+"\t    "+b_coeff
                                else: 
                                    dump_str="      "+b_atomid_n_basisid+"\t     "+b_coeff
                                #print ("when imo == p_imo", imo, dump_str)
                                imoc_valuable.append(dump_str)
                            ### Treat each imo block in sorted MOC block (6 imo block)
                            ### if index MO is changed, dump the previous imo and  initialize for the next imo
                            ### b_imo advanced already
                            #### SAVE for the previous previous imo: p_imo, p_moe
                            else:  #b_imo != p_imo:
                                #### FILTER 5 :: the imo-block will be dumped?
                                #### imo_select returns 1 for save for plot or 0 for skip
                                if vprint_log >=1 and n_check==4: print(f"line-k:{k}-{n_check}: in for-loop, when imo advances, record it");n_check+=1
                                if vprint_moc >=2: print(f"imo {p_imo}")
                                imo_select, base_coeff_dict=imo_basis(imoc_valuable, input_atoms_id, mo_type, p_imo)
                                #imo_select=f_imo_dump(imoc_valuable, atoms_id, mo_type, p_imo)
                                #### RECORD one imo if valuable 
                                if imo_select:
                                    imo_sel_ab.append(p_imo)
                                    ene_sel_ab.append(p_moe)
                                    a = QC_IMO(p_imo, p_moe, base_coeff_dict)
                                    l_QCIMO.append(a)
                                    if vprint_moc >=1:
                                        #print(base_coeff_dict)
                                        new_dict=sorted(base_coeff_dict.items(), key=lambda x:x[1], reverse=True)
                                        print(f" {p_imo} imo:: base_coeff_dict {new_dict}")
                                    #print (fl_imo_sel_ab)
                                ############# 
                                ### reset for next imo, SAVE new b_imo line Here
                                imoc_valuable=[]
                                if vprint_log >=1 and n_check==5: print(f"line-k:{k}-{n_check}: as imo advances, reset imoc_valuable[], pre_mo"); n_check+=1
                                if float(b_coeff) < 0:
                                    dump_str="      "+b_atomid_n_basisid+"\t    "+b_coeff
                                else: 
                                    dump_str="      "+b_atomid_n_basisid+"\t     "+b_coeff
                                imoc_valuable.append(dump_str)
                                p_imo = b_imo                         # p_imo advances here to present imo
                                p_moe = b_mo_ene
                        ### treat gathered last imo(6-th imo) in a block
                        ### when i_sortb == len(sorted_moc_block): last element is read, treat last b_imo in 6 imos
                        else:
                            if vprint_log >=1 and n_check==6: print(f"line-k:{k}-{n_check}: record last imo out of for-loop in sorted_moc_block");n_check+=1
                            if vprint_moc >=2: print(f"imo {p_imo}")
                            imo_select, base_coeff_dict=imo_basis(imoc_valuable, input_atoms_id, mo_type, p_imo)
                            if imo_select:
                                imo_sel_ab.append(p_imo)
                                ene_sel_ab.append(p_moe)
                                a = QC_IMO(p_imo, p_moe, base_coeff_dict)
                                l_QCIMO.append(a)
                                if vprint_moc >=1: print(f" {p_imo} imo:: base_coeff_dict {base_coeff_dict}")
                            ### last 
                            #### START NEW BLOCK :: initialize for next block
                            #imoc_valuable=[]
                            #p_imo=b_imo+1         # match imo in the next block
                            #k=0                         # k counts one block: 0 ~ nbasis+1
                            #i_block+=1
                            #moc_block=[]
                            #continue    # if k == 1+nbasis: set k=0, do not raise by next line
                            #### ? START of new IMOsave writing in advance for the next lists of imo_basis
                            #dump_str="  MO level  "+ str(b_imo) + "   "+b_mo_ene
                            #imoc_valuable.append(dump_str)
                            #if float(b_coeff) < 0:
                            #    dump_str="      "+ b_atomid_n_basisid +"\t    "+b_coeff
                            #else:
                            #    dump_str="      "+ b_atomid_n_basisid +"\t     "+b_coeff
                            #imoc_valuable.append(dump_str)
        # if flag_block == 'no': find line 
        # elif flag_block == "YES":
            # if inside block: k <= 1+nbasis
                # if k == 1+nbasis: IF LAST LINE IN MOC_BLOCK
                    # for  sorted_moc_block
                    # for-loop ends

                    ########## RESET NEW MOC BLOCK or FINISH BLOCK: initialize moc variable for new block: these cannot be used here
                    if Tag_fine == "OK":            # if you are here w imo > ilumo_cut
                        if beta == 1:
                            if tag_beta=='ON':
                                break               # if beta is done, finish
                            else:                   # if tag_beta == 'OFF': not done yet, start again
                                flag_block='OFF'    # try to find beta again
                                tag_beta='ON'
                                Tag_fine = 0
                                continue
                        else:
                            break;                  # finish here !!
                    k=0                             # k counts one block: 0 ~ nbasis+1
                    i_block+=1
                    moc_block=[]
                    continue    # if k == 1+nbasis: set k=0, do not raise by next line
                ########## CONTINUE TO NEXT LINE in the MOC block
                k += 1
                #continue
    ##print (L_imo_select, L_ene_select,L_bcdic_select, l_QCIMO, "in {}".format(whereami()))
    return L_imo_select, L_ene_select, L_bcdic_select, L_QCIMO, max_nlumo, input_atoms, input_atoms_id

###################################################### Atom symbols ###############################################
# Atom id such as C-1, C-2 are not coded yet
#FL_Atom_types=[             
#    [['Ni','d']],           # 1st type
#    [['C','p'],['O','p']]   # 2nd type
#    ]
def make_atom_basis_id(lsymbol):
    """
    make atom id and basis id
    atom id: Ni1, C1, O1, O2, H38
    basis id: s, px, dxx, f1, f2, etc
    """
    id_atom=lsymbol[0]
    atom_name=lsymbol[0]
    if re.match("\d",lsymbol[1]):
        id_atom+=lsymbol[1]        # id_atom == Ni1, Ni2, H30, H31, etc
        l_basis=lsymbol[2:]
    else:
        #id_atom+='1'               # id_atom == Ni C1 C O1 O2
        l_basis=lsymbol[1:]
    id_basis=l_basis[0]
    if len(l_basis) == 2:
        id_basis+=l_basis[1]
    return atom_name, id_atom, id_basis


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
    print(latom_name, latom_index, "in", whereami())
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

    print(f"all atoms in {whereami()}(): {all_atom_count}")
    return all_atom_count

def Cal_Ncore(dict):
    ncore = 0
    for atom in dict:
        nc = Atom_Core_631g[Atom_Table.index(atom)]
        ncore += nc * dict[atom]
    print(f"all the number of core basis in {whereami()}(): {ncore}")
    return ncore        


####################### imo_basis #########################################################    
def imo_basis(imoc_list, l_atoms, filter_tag, mo_id):
    """
    imoc_list:: "atomid_and_base, moceff"
    decide whether 1 level of dump block of 'imoc_list' as for the given mo_id is saved for draw
    imoc_list is saved for the atoms listed in atom_list
    """
    
    #### change element which corresponds to l_atoms ordering, to check its existence with high MO coefficients
    n_latoms = len(l_atoms)
    atom_tag=np.zeros(n_latoms)
    base_coeff = {}
    mcoeff = 0.0
    #print l_atoms, "in", whereami()
    if not imoc_list:
        return 0, base_coeff
    elif l_atoms == None:
        return 1, base_coeff
    else:
        k=0
        ### Obtain atom_tag: which atom is included in imo one block [0,0,0] => [1,0,1] 1st and 3rd atoms exist
        for atom in l_atoms:
            flag="OFF"
            j=0
            for moc_line in imoc_list:
                #### skip first line of "MO level  mo-id  mo-ene"
                #print moc_line, "in", whereami()
                if j==0:
                    j+=1
                    continue
                #### in C-1-py 0.27933, we can find atom species in capital C,
                #### py, s, dx, f etc is lowercase
                #### do not use search for N as for "N Ni"
                fields=moc_line.split()
                fields=[x for x in fields if x]
                t_atom_basis=re.split("-", fields[0])
                atom_id=t_atom_basis[0]
                #if re.search(atom, moc_line):
                #print atom, atom_name
                ### atom existence is checked in atom_tag
                if atom == atom_id:
                    flag="ON"
                    #print atom, " in ", moc_line
                    break
            if flag=="ON":
                atom_tag[l_atoms.index(atom)]=1
        
        #   for SUB, atom_tag == [1,0,0,1] etc: 1st atom should exist
        #print atom_tag, "in", mo_id, "with", filter_tag, "in", whereami()
        if filter_tag=="ONE" :
            if atom_tag[0]==1:
                base_coeff = imo_dic_basis_coeff(imoc_list)
                selection= 1
            else: 
                selection= 0
        elif re.search("SUB", filter_tag):
            if atom_tag[0]==1: # Ni is neccessary:: this fails of CO2 MO alone in complex
                base_coeff = imo_dic_basis_coeff(imoc_list)
                selection= 1
            else:
                selection= 0
        elif re.search("SEL", filter_tag):
            if 1 in atom_tag:
                base_coeff = imo_dic_basis_coeff(imoc_list)
                selection= 1
            else:
                selection= 0
        elif filter_tag=="ALL":             
                base_coeff = imo_dic_basis_coeff(imoc_list)
                selection= 1
        else:
            print ("Error:: No atom filtering type in", whereami() )
            exit(55)
    #print base_coeff
    #new_dic = trim_coeff(base_coeff, n_latoms*Nbasis_show)
    new_dic = trim_coeff(base_coeff, Nbasis_show)
    if vprint_moc >=2: 
        if vprint_loop==1: print (f"mo_id {mo_id} in {whereami()}()")
        print("trimmed dict:: ", new_dic)
    return selection, new_dic

def trim_coeff(dic, n):
    """ trim dict to leave some maximum values """
    if vprint_moc >= 2:
        print (whereami())
        #lprint_sorted(dic)
    if len(dic) > n:
        a = heapq.nlargest(n, dic.items(), lambda i: i[1])
        if vprint_moc >= 2: lprint(a)    
        new_dic = dict(a)
        return new_dic
    else:
        return dic

from decimal import *
getcontext().prec = 5
def imo_dic_basis_coeff(imoc_list):
    """
    imoc_list :: atomid_and_basis, mocoeff
    make basis dictionary
    Return format:: {'P2-px': 0.27074, 'P1-px': 0.27057}
    in case there are two same bases:: arithmatically add at the moment
    """
    i=0
    bcoeff={}
    for line in imoc_list:
        #if i == 0:
        #    if vprint_moc >= 1: print (f"imoc list:{line}:{whereami()}()")
        #else:
        if vprint_moc >= 3: print (f"imoc list:{line}:{whereami()}()")
        base_coeff = line.split()
        coeff = abs(float(base_coeff[1]))
        base = base_coeff[0]
        if not base in bcoeff.keys():
            bcoeff[base] = coeff
        else:
            ### arithemetic sum
            tmp = Decimal(bcoeff[base]) + Decimal(coeff)
            ### sum of square of coefficient
            tmp = (Decimal(bcoeff[base]*bcoeff[base]) + Decimal(coeff*coeff)).sqrt()
            bcoeff[base] = float(tmp)
        i+=1
    if i==1:
        pass
    #print(bcoeff)
    #sorted_bcoeff=sorted(bcoeff.items(), key=lambda x:float(x[1]))
    if vprint_moc >= 2:
        print (f"bases_coeff dictionary: {bcoeff} in {whereami()}()")
    return bcoeff




###################################### BLOCK 3 :: NBO Bonding extraction #######################
def analyze_nbo(nfile, fid):
    return 0
    """
    if nfile==3 and fid == 1:   #### As for complex molecule
        if flag_nbo == "NO" and tag_nbo == "YES":
            if re.search(KW_nbo, line):
                flag_nbo = "YES"
                print ("#############  NBO Analysis  ##########")
                continue
        #### found NBO and the next line is "---------------------------------"                
        elif flag_nbo == "YES" and tag_nbo == "YES":
            if flag_bd=="NO":
                #### Only Ni is checked
                #if re.search(KW_nbo_BD, line) and re.search('Ni', line):
                if re.search(KW_nbo_BD, line) and (re.search('Ni', line) or re.search('O', line)):
                    flag_bd="YES"
                    bd_line=re.split("[\s\(\)-]+", line)
                    bd_line=[x for x in bd_line if x]
                    print (bd_line[0], bd_line[1], bd_line[2], bd_line[3], bd_line[4], bd_line[5], bd_line[6], bd_line[7] )
                    nb_line.append(nbo_BD_nline_atoms[nbo_atoms.index(bd_line[4])])
                    nb_line.append(nbo_BD_nline_atoms[nbo_atoms.index(bd_line[6])])
                    #print nb_line
                    ibd += 1
                elif re.search(KW_nbo_end, line):
                    tag_nbo = "Done"
                else:
                    pass
                continue
            ### in the NBO-BD block, write
            else:
                if ibd == 1:
                    #print line
                    a=line[15:24]
                    b=line[27:33]
                    c=line[34:36]
                    d=line[40:41]
                    e=line[41:50]
                    f=line[50:51]
                    g=line[56:65]
                    if len(line) > 70:
                        h=line[65:66]
                        hh=line[71:80]
                        print ("\t", a, b, c, d, e, f, g, h, hh)
                    else:
                        print ("\t", a, b, c, d, e, f, g)
                    #print line
                    ibd+=1
                    continue
                if ibd < nb_line[ibd_atom]+1:
                    ibd+=1
                    continue
                else:
                    if ibd_atom==0:
                        ibd=1       #### to read second atom part
                        ibd_atom+=1
                        continue
                    else:
                        ibd=0
                        ibd_atom=0
                        nb_line=[]
                        flag_bd="NO"
                        continue
    """                        
