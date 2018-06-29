#!/usr/bin/python
### made by J. Park
### 2016. 10. 25 developed for radical
### 2016. 10. 28 developed for line
### 2016. 11. 14 developed for homo-reference
### 2016. 11. 15 developed for nbo BD
### 2016. 11. 21 developed for MO Coefficients analysis
### 2016. 11. 24 developing for nbo charge
### 2016. 12. 13 developed for MO linkage using hash
### 2017. 01. 24 developing for 6-31G(d)
### 2017. 02. 06 developing for ini.py
### 2017. 02. 14 up to 5 files
### 2018. 05. 14 atom two group for link_split
### git test
### argparse

import re
import os
import argparse
import _math
import numpy as np
from common import *
from mplt_qcdraw import *
from mp_mo_dump import *
from mp_file_anal import *
from mplt_mo_ini import *
from qmo import *

################### MAIN LOOP ####################
def analyze_files(files,atomlists,atom_group,motypes, base_type,L_link, link_types, atom_groups): 
    # start analyze_files function 
    print files, "MO filter[atom|type]", atomlists, motypes,"L_link=", L_link, "link_type=", link_types
   
    """
        FL_:: all-files variable [[[file1-a],[file1-b]],[[file2-a],[file2-b]],...
              w. a,b dependency,  last dimension is 2
        Fl_:: all-files variable [[file1],[file2],...
        fL_:: for single file variable with a, b dependency
        fl_:: for single file variable 
        FILTER 0: display for selected atoms in case -a
        FILTER 1: core orbital is excluded
        FILTER 2: cut by MO coefficient criterion
        FILTER 3: cut by max_nlumo
        FILTER 4: 
    """
    Fl_qchem_outf = files[:]
    nfile=len(Fl_qchem_outf)
    #### FILTER MO for selected atoms
    if atomlists:                   # go to MOC BLOCK 2
        L_mo_coeff = True
        MOcoeff_crit=0.15           # 0.15, 0 for all the levels
    else:
        L_mo_coeff = False
        max_lumo=Nmax_virtual
    

    ### MO block
    FL_norbital=[]
    FL_homo_ene=[]      # it can have beta, save homo energy [ [f1 a, f1 b],[f2 a, f2 b]... ]
    FL_homo_id=[]       # ene and MO id go parallel
    Fl_homo_ene=[]
    Fl_beta_tag=[]      # 0 for restricted, 1 for beta spin of  unrestricted calculation 


    ### MOC block
    Fl_MO_types=[]      # [ONE(f1), 1sub(f2)...
    Fl_MO_atoms=[]      # just different atoms: [f1: [Ni, P], f2: [C, O], ..
    Fl_MO_atoms_id=[]   # all indices
    Fl_nlumo=[]         # [nlumo-f1, nlumo-f2...    same as a|b
    Fl_QCIMO=[]         # class list
    FL_MOene_list=[]      # [ [e1, e2... ], [e1, e2...] [...
    FL_MOene_select=[]
    FL_MObcdic_select=[]  # base-coefficient selected
    FL_abMO_dic_list=[]
    FL_MOid_select=[]
    #### depending on nbasis, num of nbo BD line is different
    nbo_atoms=['Ni', 'P', 'C','O']    # atom for 1st fragment, 2nd fragment
    nbasis=[22, 8, 9, 9]
    nbo_BD_nline_atoms=[]
    for nb in nbasis:
        nb_line5=int(nb/5)+1
        nbo_BD_nline_atoms.append(nb_line5)

    KW_nao_occ="NATURAL POPULATIONS"
    KW_nec="Natrual Electron"           # This has information of num of e in NAO
    KW_nbo="Bond orbital"
    KW_nbo_BD="BD"
    KW_nbo_end="NBO state"

    AngSym=('1s', '2s', '2px', '2py', '2pz')

    lcmo_mol=1      
    lcmo_nlnh=[1,1]
    lcmoA=[]
    lcmoB=[]
    lcmoC=[]

    ### LINK block
    FL_sMO_dic=[]
    FL_link_dic=[]
    nlink =  NLINK
    """
    Constants in module
    readline file and extract (1) energies, (2) coefficients, (3) nbos
    """
    fid=0
    s='-'       # atom-id linker
    ### Analyze each file
    #####################
    for file in Fl_qchem_outf:
        # tags to get MO energies
        tag_MO="YES"    # YES, DONE, FINISH MO ene
        flag_MO="OFF"
        f_homo_i_line=[]        # one or two(beta_tag) values of index of homo line
        f_line_ene_a=[]     # for job=optimization, whenever find new keyword, empty
        f_line_ene_b=[]
        fL_line_ene=[]
        fL_homo_id=[]
        nbasis=0

        All_atoms=[]
        
        # tags to get MO coefficients
        if L_mo_coeff:
            tag_moc="YES"
        else:
            tag_moc="NO"
        flag_moc="OFF"
        flag_moc_beta="OFF"
        temp_atom_name='Z'    # not atom id such as Ni, Ni-1, Ni-2...
        imoc=0
        imoc_block=0
        imo_tmp=0       # mo-id starts from 1 in MO Coefficient Block
        fl_QCIMO={}     # class list

        # tags for nbo
        tag_nbo="YES"
        flag_nbo="NO"
        flag_bd="NO"

        fl_ene_sel_ab=[]
        fl_imo_sel_ab=[]
        fl_bcdic_sel_ab=[]
        fl_ene_select=[]
        fl_imo_select=[]
        fl_bcdic_select=[]
        f_ene_dic={}
        f_ene_dic_beta={}
        f_dic_ab_list=[]
        f_norbital_ab=[]
        disp_MOids=[]
        ibd=0
        ibd_atom=0
        nb_line=[]
        mo_ene_tmp=-100
        imoc_valuable=[]

        # read one file for analysis & write tmp_file.txt
        #inf=open(file, 'r')
        lfile = re.split('\.', file)
        tmp0 = 'tmp0-' +lfile[0]+'.log'
        tmp1 = 'tmp1-' +lfile[0]+'.log'
        if fid == 2:
            fname_link = lfile[0]   # save to write link pair

        i=0
        ib_eline=0     # beta ene line index
        ### Read 1 File
        ###############
        #while 1:
        with open(file, 'r') as f:
            for line in f:      #=f.readline()
                i += 1
                if not line: break         # end of while-loop, one file
                ### Find nbasis for MO Coefficient in advance
                if L_mo_coeff and nbasis==0:
                    if re.search(KW_nbasis, line) and re.search("shells", line):
                        field=re.split("\s+", line)
                        field=[x for x in field if x]
                        nbasis=float(field[5])
                        #print "found num of basis: ", nbasis
                ### BLOCK 1 :: MO energy block
                ### if-MO block
                if flag_MO == "OFF" and tag_MO == "YES":
                    ### Find MO block & initialize (first "occupied" is for Alpha spin)
                    if re.search(KW_MOene_occ, line):
                        print "### ", Fl_qchem_outf[fid], "#########################################"
                        print "###  BLOCK 1 :: QChem MO energy"
                        flag_MO="ON"
                        Fl_beta_tag.append(0)   # beta option is initiated
                        ### initiate all the variables used in this block
                        i_ene_line=0            # energy line index starts from 0
                        tag_beta="NO"
                        flag_vline="NO"         # now occupied orbital
                        continue
                #### inside MO Block
                elif flag_MO == "ON" and tag_MO == "YES":
                    ### if options except energies, find keyword in reverse order
                    ### 1. find end block line: if arrived the END line of MO energy
                    if re.search(KW_MOene_end, line):
                        tag_MO="Done"   # escape MO block
                        continue
                    #### 2. skip non-ene line (symmetry lines: watch out B|Beta)
                    elif re.search("A", line) or (re.search("B", line) and not re.search("Beta", line)):
                        continue
                    #### 3. find Virtual line
                    elif re.search(KW_MOene_vir, line):
                        flag_vline="YES"    # Virtual orbital energies
                        #   for H atom (one electron - no beta occ)
                        if Fl_beta_tag[fid] == 1:
                            if not f_line_ene_b:
                                print f_line_ene_b, "is false"
                                iline_homo=-1
                        f_homo_i_line.append(iline_homo) # twiced if beta_tag, iline_homo is assigend below
                        print "f_homo_i_line ", f_homo_i_line
                        continue
                    #### in case of Beta, there is blank line, how to treat BLANK
                    #### 4. find blank for beta
                    elif not line.strip():
                        flag_vline="NO"     # virtual energy ended and prepare for beta
                        continue
                    #### 5. find BETA
                    elif re.search(KW_MOene_beta, line):
                        tag_beta="YES"
                        Fl_beta_tag[fid]=1  # BETA, unrestricted calculation is corrected here
                        i_ene_line=0        # reset energy line index, it starts from 0
                        continue
                    ### 6. find "Occupied" again for beta
                    elif re.search(KW_MOene_occ, line):
                        continue
                    #### save energy line if skipped all the previous if condition
                    else:
                        if tag_beta == "YES":
                            ### if beta is detected but is the same as alpha, igmore for bent CO2 in 1-PP-B
                            if ib_eline == 0:
                                print "check beta line", line.rstrip, f_line_ene_a[ib_eline]
                                if line.rstrip() == f_line_ene_a[ib_eline]:
                                    print "Beta spin is detected but same as Alpha spin"
                                    Fl_beta_tag[fid]=0
                                    tag_MO="Done"
                                    continue
                            f_line_ene_b.append(line.rstrip())
                            if flag_vline != "YES":
                                iline_homo=i_ene_line
                            ib_eline+=1
                        else:
                            f_line_ene_a.append(line.rstrip())
                            if flag_vline != "YES":
                                iline_homo=i_ene_line
                        i_ene_line += 1
                        continue
                ### if-MO block
                ### Wrap up MO Block
                elif tag_MO == 'Done':
                    """
                        f_line_ene_a=[l1, l2...]
                        fL_line_ene =[[Alpha:l1, l2,...],[Beta:l1, l2,...]]
                        FL_homo_ene =[[f1:[A-homo, B-homo],[f2:...
                        Fl_homo_ene =[f1-homo, f2-homo,...]
                        FL_homo_id  =...
                    """
                    tag_MO = 'Finish'
                    print tag_MO
                    fL_line_ene.append(f_line_ene_a)
                    if Fl_beta_tag[fid] == 1:
                        fL_line_ene.append(f_line_ene_b)
                    fL_ene_ab,f_homo_ab,f_homo_id_ab=f_extract_ene(f_homo_i_line, fL_line_ene)
                    FL_homo_ene.append(f_homo_ab)
                    print FL_homo_ene, f_homo_id_ab
                    if Fl_beta_tag[fid] == 1:
                        Fl_homo_ene.append( max(FL_homo_ene[fid][0], FL_homo_ene[fid][1]) )
                    else:
                        Fl_homo_ene.append( FL_homo_ene[fid][0])
                    FL_homo_id.append(f_homo_id_ab)
                    print FL_homo_id
                    FL_MOene_list.append(fL_ene_ab)
                    
                    f_norbital_ab.append(len(fL_ene_ab[0]))
                    if Fl_beta_tag[fid] == 1:
                        f_norbital_ab.append(len(fL_ene_ab[1]))
                    FL_norbital.append(f_norbital_ab)
                    continue
                
        inf=open(file, 'r')
        while 1:    #(realline file)
        #with open(file, 'r') as inf:
            line=inf.readline()
            
            ### BLOCK 2 :: MO Coeff ########################################
            ### ene & index are obtained here
            ### find key_word_MOcoeff and initialize
            if flag_moc == "OFF" and tag_moc=="YES" :
                if re.search(KW_MOcoeff, line) :    #common for alpha & beta
                    flag_moc="ON"
                    moc_block=[]
                    #### constants from module mplot_ini
                    max_nlumo = Nmax_virtual
                    Fl_nlumo.append(max_nlumo)
                    #### obtain atoms for MO and mol types in QChem outfile
                    if atomlists:
                        fl_mo_atoms_id=atomlists[fid].split()
                        atom_style='name'
                        for atomid in fl_mo_atoms_id:
                            if re.search('\d', atomid):
                                atom_style='index'
                        if atom_style=='name':
                            fl_mo_atoms=fl_mo_atoms_id[:]
                            Fl_MO_atoms.append(fl_mo_atoms)
                            Fl_MO_atoms_id.append(fl_mo_atoms)
                        ### USE index
                        elif atom_style=='index':
                            latom_name, latom_index,atom_species= atom_decompose(fl_mo_atoms_id)
                            fl_mo_atoms=latom_name[:]
                            #fl_mo_atoms=atom_species[:]
                            Fl_MO_atoms.append(fl_mo_atoms)
                            Fl_MO_atoms_id.append(fl_mo_atoms_id)
                            #Fl_MO_atoms.append(atom_species)
                    if motypes:
                        f_mo_type=motypes[fid]
                    else:
                        f_mo_type=get_atom2motype(Fl_MO_atoms[fid])         
                    Fl_MO_types.append(f_mo_type)             #Fl_MO_type[:]
                    print fl_mo_atoms_id, Fl_MO_atoms_id, f_mo_type, Fl_MO_types, "in", whereami()
                    #Fl_MO_atoms, Fl_MO_type = get_mo_labels(files, atoms, i_mo_types)
                    #Fl_nlumo, MO_link_hl_id = get_nlumo_linkid(nfile)
                    #print fl_mo_atoms, f_mo_type, Fl_nlumo, MO_link_hl_id
                    ### use common homo index for a|b
                    fL_homo_id=FL_homo_id[fid]      # list of one|two elements
                    f_homo_a = FL_homo_id[fid][0] # Alpha has always large value
                    f_lumo_id_cut=f_homo_a+max_nlumo
                    #lumo_id_cut=f_homo_id+Fl_nlumo[fid]
                    if f_norbital_ab < f_lumo_id_cut:
                        f_lumo_id_cut = f_norbital_ab
                    print "MO cut id: ", fL_homo_id, f_lumo_id_cut
                    if flag_moc_beta == "ON":
                        f_homo_id=FL_homo_id[fid][1]
                        print "###################  BLOCK 2 :: MO Coeff beta starts "
                        imoc=0
                        imoc_block=0
                        imo_tmp=0            # mo-id starts from 1 in MO Coefficient Block
                        mo_ene_tmp=-100
                        imoc_valuable=[]
                        fl_imo_sel_ab=[]
                        fl_ene_sel_ab=[]
                    else:
                        print "###################  BLOCK 2 :: MO Coeff starts "
                    print Fl_qchem_outf[fid], ":: MO Levels"
                    continue
            ### Inside MOC PART
            elif flag_moc == "ON" and tag_moc=="YES":
                #### if end of MO Coefficients
                if not line.strip():
                    if Fl_beta_tag[fid]==0:
                        print "End of MOC of alpha"
                    continue
                #### First MOC-block line: save MO id
                if imoc==0:
                    imo=line.split()
                    if imoc_block == 0:
                        Nene_mocline=len(imo)
                    if V_print_moc >=2: print "imoc==0; n_ele_col=", Nene_mocline, line
                    imoc+=1
                    continue
                #### save energy
                elif imoc==1:
                    imo_ene=line.split()
                    imo_ene.pop(0)      # first column is "eigen values"
                    if V_print_moc >= 2: print imo_ene
                    imoc+=1
                    continue
                ###### analyze ONE MOC line for one basis in MOC_BLOCK
                elif imoc <= 1+nbasis:      # nbasis can be got in the first part before MO block
                    # moc_line_analysis(line, Nene_mocline, temp_atom_name)
                    lmoc_line=line.split()
                    if V_print_moc >=2: print lmoc_line
                    ncol = len(lmoc_line)
                    ind = ncol - Nene_mocline
                    lcoeff = lmoc_line[ind:]
                    l_symbol = lmoc_line[0:ind]
                    ibasis = l_symbol.pop(0)
                    #print l_symbol, whereami()
                    ### make atom id
                    atom_name, id_atom, id_basis = make_atom_basis_id(l_symbol)
                    if not id_atom in All_atoms:
                        All_atoms.append(id_atom)
                    #if V_print_moc >= 1: print ibasis, l_symbol, lcoeff, atom_name, id_atom
                    ####  Each MOC line is decomposed to maximum Nene_mocline==6 basis-ene-coeff data
                    for j in range(Nene_mocline):
                        mo_symbol=[]
                        moc_tmp=[]
                        #### FILTER 1: 6 data of MOC line filtered by MOcoeff_crit
                        if abs(float(lcoeff[j])) > MOcoeff_crit:
                            mo_symbol.extend((id_atom, id_basis))
                            moc_tmp.extend((mo_symbol,imo[j],imo_ene[j],lcoeff[j]))
                            moc_block.append(moc_tmp)   # moc_block contains required data type for 1 MOC block ( Nene_mocline*nbasis maximum)
                            #print(np.matrix(moc_block))

                    ####### Treat one MOC block (6 mo_id) and INITIALIZE variables again
                    if imoc == 1+nbasis:
                        ### for first MOC block, calculate core
                        if imoc_block == 0:
                            all_atom_count = All_atom_count(All_atoms)
                            ncore = Cal_Ncore(all_atom_count)

                        ### initialize moc variable for new block
                        imoc=0
                        imoc_block+=1
                        #### treat block of MOC here
                        #### new style:           [['Ni4', 'd4'], '100', '-0.016', '0.19560']
                        #### FORMAT of moc_block: [[atom_id,basis_name], mo_id, mo_ene, mo_coeff
                        #### SORTED by  mo_id
                        sorted_moc_block=sorted(moc_block, key=lambda x:float(x[1]))
                        if V_print_moc >=2: 
                            for  x in sorted_moc_block:
                                print x
                        #### scan sorted MOC block: by one line - format
                        if V_print_moc >= 2: print "Loop to sorted MOC line"
                        for eff_coeff in sorted_moc_block: # to the end
                            ### FILTER 2  skip core orbital based on atom species in lmoc_line[1]: follow QChem basis
                            if int(eff_coeff[1]) <= ncore:            # NCORE is treated here
                                continue
                            ### FILTER 3: cut lines using FL_nlumo[0] 
                            elif int(eff_coeff[1]) > f_lumo_id_cut:
                                #### dump lastly MO Coefficient before break
                                if V_print_log >=1: print "last: ", fl_imo_sel_ab
                                imo_select,bas_co_dic=f_imo_basis(imoc_valuable, Fl_MO_atoms_id[fid], Fl_MO_types[fid],imo_tmp)
                                #imo_select=f_imo_dump(imoc_valuable, Fl_MO_atoms[fid], Fl_MO_types[fid],imo_tmp)
                                if imo_select:
                                    fl_imo_sel_ab.append(imo_tmp)
                                    fl_ene_sel_ab.append(mo_ene_tmp)
                                    a = QC_imo(mo_ene_tmp, base_coeff_dict)
                                    fl_bcdic_sel_ab.append(a)
                                    fl_QCIMO[imo_tmp]=a
                                #Fl_QCIMO.append(fl_QCIMO)
                                #### save imo & ene for alpha or beta
                                fl_imo_select.append(fl_imo_sel_ab)
                                fl_ene_select.append(fl_ene_sel_ab)
                                fl_bcdic_select.append(fl_ene_sel_ab)
                                ### Finalize MOC
                                #print fl_imo_select[0]
                                if Fl_beta_tag[fid]==0:
                                    print "Cut by nlumo"
                                    tag_moc="Done"
                                elif Fl_beta_tag[fid]==1 and flag_moc_beta=="ON":
                                    tag_moc="Done"
                                #### run beta again
                                elif Fl_beta_tag[fid]==1 and flag_moc_beta=="OFF":
                                    flag_moc_beta="ON"
                                    flag_moc="OFF"       # find KW for beta
                                #### FINISH MO Coeff selection
                                break
                            if tag_moc=="Done":
                                continue
                            ######## DUMP for selected atom line in selected 6 mo-id's ##########
                            ####### FILTER 4       :: cut by atom species and ID here
                            ####    f_imo_dump()   :: decision for dumping the imo-block
                            ####    eff_coeff== [['Ni4', 'd4'], '100', '-0.016', '0.19560']
                            #### 'pass' writes the line in dump string      of str_dump
                            #### 'continue' skips the line before writing   in str_dump
                            if Fl_MO_types[fid] == 'ALL':
                                pass
                            elif Fl_MO_types[fid] == "ONE":
                                if eff_coeff[0][0] == Fl_MO_atoms_id[fid][0]:
                                    if V_print_filter >= 2: print eff_coeff[0][0], Fl_MO_atoms_id[fid][0], "in filter 4"
                                    pass
                                else:
                                    continue
                            #### 'ALL' and '1sub' works the same but different in f_imo_dump()
                            #### more detail screening difference between 'ALL' and '1sub', refer to f_imo_dump
                            #### for more options for atom id w. index
                            elif (re.search("sel",Fl_MO_types[fid],re.IGNORECASE) or re.search("sub",Fl_MO_types[fid], re.IGNORECASE)):
                                if eff_coeff[0][0] in Fl_MO_atoms_id[fid]:
                                    pass
                                else:
                                    continue 
                            else:
                                print "Error:: there is no atom type in", whereami() 
                                exit(11)
                            #print eff_coeff
                            imo=int(eff_coeff[1])
                            mo_ene=eff_coeff[2]
                            #### if index MO is the same, save for write
                            if imo == imo_tmp:
                                if float(eff_coeff[3]) < 0:
                                    dump_str="      "+s.join(eff_coeff[0])+"\t    "+eff_coeff[3]
                                else: 
                                    dump_str="      "+s.join(eff_coeff[0])+"\t     "+eff_coeff[3]
                                #print imo, dump_str
                                imoc_valuable.append(dump_str)
                            ### Treat each imo block in sorted MOC block (6 imo block)
                            ### if index MO is changed, dump the previous imo and  initialize for the next imo
                            else:
                                #### dump for one-imo, the previous saved moc's
                                #### FILTER 5 :: the imo-block will be dumped?
                                #### imo_select returns 1 for save for plot or 0 for skip
                                #imo_select,mbasis_dic=f_imo_basis(imoc_valuable, Fl_MO_atoms[fid], Fl_MO_types[fid], imo_tmp)
                                imo_select, base_coeff_dict=f_imo_basis(imoc_valuable, Fl_MO_atoms_id[fid], Fl_MO_types[fid], imo_tmp)
                                #imo_select=f_imo_dump(imoc_valuable, Fl_MO_atoms_id[fid], Fl_MO_types[fid], imo_tmp)
                                #### RECORD if valuable 
                                if imo_select:
                                    fl_imo_sel_ab.append(imo_tmp)
                                    fl_ene_sel_ab.append(mo_ene_tmp)
                                    a = QC_imo(mo_ene_tmp, base_coeff_dict)
                                    fl_bcdic_sel_ab.append(a)
                                    fl_QCIMO[imo_tmp]=a
                                    #print fl_imo_sel_ab
                                #### START of NEW IMO :: initialize for next imoc block
                                imoc_valuable=[]
                                imo_tmp=imo
                                mo_ene_tmp=float(mo_ene)
                                #### START of new IMOsave writing in advance for the next lists of imo_basis
                                dump_str="  MO level  "+ str(imo) + "   "+mo_ene
                                imoc_valuable.append(dump_str)
                                atom_basis=s.join(eff_coeff[0])
                                if float(eff_coeff[3]) < 0:
                                    dump_str="      "+ atom_basis +"\t    "+eff_coeff[3]
                                else:
                                    dump_str="      "+ atom_basis +"\t     "+eff_coeff[3]
                                imoc_valuable.append(dump_str)
                                ### for link: basis control - basis can be repeated
                                #a = QC_imo(imo_tmp, atom_basis, mo_ene_tmp)

                        #### end of MO coeff one block
                        moc_block=[]
                        continue
                    #### have read one block of nbasis lines == MOC block
                imoc+=1
                continue
        # while
            ##### BLOCK 3 :: NBO BonDing extraction #################
            if L_NBO_analysis == 1:
                if nfile==3 and fid == 1:   #### As for complex molecule
                    if flag_nbo == "NO" and tag_nbo == "YES":
                        if re.search(KW_nbo, line):
                            flag_nbo = "YES"
                            print "#############  NBO Analysis  ##########"
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
                                print bd_line[0], bd_line[1], bd_line[2], bd_line[3], bd_line[4], bd_line[5], bd_line[6], bd_line[7] 
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
                                    print "\t", a, b, c, d, e, f, g, h, hh
                                else:
                                    print "\t", a, b, c, d, e, f, g
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
# func
    # for file
        # while
        # end of read 1 file
        inf.close() 
        ######## FILTER 5: retailing plot list in each file
        #### cut lower energy limit for drawing in figure
        #### make dictionary for each file
        #### change into all levels 
        print "#### End of file (close file):: Filter 5 #######"
        if not atomlists :
            not_digit=[]
            i_lumo_cut= min(f_norbital_ab[0], FL_homo_id[fid][0]+Nmax_virtual)
            i_homo_cut= max(1,FL_homo_id[fid][0]-Nmax_occupy)
            print "draw all level :", i_homo_cut, i_lumo_cut
            if V_print == 1: print fL_ene_ab
            for ab_line1d in fL_ene_ab:    # or FL_homo_id[fid]
                tmp_list=[]
                tmp_elist=[]
                j=1
                for ene in ab_line1d:
                    if re.search("\*", ene):
                        pass
                    # select using index                        
                    elif i_homo_cut <= j and j <= i_lumo_cut:
                        #print ene
                        tmp_elist.append(float(ene))
                        tmp_list.append(j)
                    j+=1
            
                fl_ene_select.append(tmp_elist[:])
                fl_imo_select.append(tmp_list[:])

        if V_print >= 3:
            print "After imo selection, make dictionary: ", fl_imo_select
            print fl_ene_select
        print "<Plot MO levels of ", Fl_qchem_outf[fid], ">"
        print "   MO_level_ID Energy"
        for imo, moe in zip(fl_imo_select[0],fl_ene_select[0]) :
            #i_imo=imo
            #f_mo=moe
            print "%11d %8.3f" % (imo, moe)
            f_ene_dic[imo]=moe
        f_dic_ab_list.append(f_ene_dic)
        if Fl_beta_tag[fid] == 1:
            #### change key for alpha, beta mo-id
            print fl_imo_select[1]
            print fl_ene_select[1], len(fl_ene_select[1])
            for imo, moe in zip(fl_imo_select[1],fl_ene_select[1]) :
                i_imo=imo
                f_mo=moe
                print i_imo, moe
                f_ene_dic_beta[i_imo]=f_mo
            f_dic_ab_list.append(f_ene_dic_beta)
        #print f_dic_ab_list            
        #### code for beta is developed up to here
        #print f_ene_dic
        FL_MOid_select.append(fl_imo_select[0])
        FL_MOene_select.append(fl_ene_select[0])
        FL_MObcdic_select.append(fl_bcdic_select[0])
        FL_sMO_dic.append(f_ene_dic)
        FL_abMO_dic_list.append(f_dic_ab_list)
        fid += 1
    ######################################################
    ###### end of reading all  the files
    print "############################################"
    print "############### Result #####################"
    print "       Filename     [HOMO_energy] [ID] beta_tag motype linktype"
    for i in range(nfile): 
        print "     %-15s" % Fl_qchem_outf[i], 
        print FL_homo_ene[i],
        print FL_homo_id[i],
        print Fl_beta_tag[i],
        print motypes[i] 

    ################# Block 5 ::PLOT using MPL ###########
    #str="--------"
    """
    MO_link_id=(94,99)
                        for left-side link      for right-side link
                     1st-mol 2nd-mol            2nd-mol 3rd-mol
    MOlink_hl_id=[[[ 0       ,0  ],[0 ,   1]], [[0,       1],[    1,     1]]]
    We need          y-val  y-val y-val y-val 
    """
    print "##### Block5: Final MO energies for mpl plot:: "

    MOlink_hl_id=[[[ 0       ,0  ],[1 ,   1]]]

    """
    Link is only working for non-degenerate MOs
    """
    MO_link_hl_ene=[]

    ######## Make MO-link pair
    ######## how to obtain link pair

    #### To use mpl-plot, it needs energies for y-axis, not id
    if L_link:
        #if 'MOlink_hl_id' not in locals():
        #    print "Error:: No link pair input in", os.path_basename(__file__)
        #    exit(90)
        print "############ Link Plot #####################"
        print len(Fl_QCIMO)
        for fl_class in Fl_QCIMO:
            for qc_class in fl_class:
                attrs = vars(qc_class)
                print ', '.join("%s: %s" % item for item in attrs.items()) 
        ### LOOP for nfile-1 links                
        for i in range(nfile-1): 
            #iside_link=0     # 0 for between 0th and 1st files
            ihomo_f0=FL_homo_id[i][0]       # f0 for left file, f1 for right,
            ihomo_f1=FL_homo_id[i+1][0]     # due to alpha, beta
            #### here obtain link pair
            print "Link ID in %d-th interval:: " % (i+1)  
            #### find MO id using homo-id
            #### between MO 0 and MO 2
            #### obtain link pair here
            #f_link
            link_y=[]       # make link pair
            
            if link_types[i] == 'flow':
                fl_molink_hl_id=get_link(Fl_QCIMO[i],ihomo_f0,Fl_QCIMO[i+1],ihomo_f1,'flow', nlink, atom_groups)
            ### left mo has lower number
            elif link_types[i] == 'lsplit':
                fl_molink_hl_id=get_link(Fl_QCIMO[i],ihomo_f0,Fl_QCIMO[i+1],ihomo_f1,'lsplit', nlink, atom_groups)
            elif link_types[i] == 'rsplit':
                fl_molink_hl_id=get_link(Fl_QCIMO[i],ihomo_f0,Fl_QCIMO[i+1],ihomo_f1,'split', nlink, atom_groups)
            ### write MO Link pair to repair and reuse it
            #with open("link_id.dat", "w") as f:
            #    f.write fl_molink_hl_id

            for l_pair in fl_molink_hl_id:
                print l_pair, 
                key_f0=ihomo_f0+l_pair[0]
                key_f1=ihomo_f1+l_pair[1]
                if key_f0 not in FL_sMO_dic[i]:
                    print "Link dictionary error:: exit(44)"
                    print "use option L=0"
                    exit(44)
                if key_f1 not in FL_sMO_dic[i+1]:
                    print key_f2, FL_sMO_dic[i+1]
                    print "Link dictionary error:: exit(45)"
                    print "use option L=0"
                    exit(45)
                l_ene1=FL_sMO_dic[i][key_f0]
                l_ene2=FL_sMO_dic[i+1][key_f1]
                link_ene=[l_ene1, l_ene2]
                link_y.append(link_ene)
            print " "
            MO_link_hl_ene.append(link_y)
            #print link_y
        #### 3D list, draw with x-coord's
        if V_print_link >= 0:
            print "Link energy:: "
            print MO_link_hl_ene
    """
        type ab_draw(int, int, int, int, hash(int, float), float)
            sub function: 
    """

    fid=0
    #### There might be alpha, beta two dictionaries for one file
    for f_abMO_dic in FL_abMO_dic_list:
        #### draw one or two dictionary
        beta_tag=Fl_beta_tag[fid]
        if beta_tag==0:
            # number of files, file id, beta?, alpha or beta, dic[a or b], homo(a or b)
            mplt_ab_draw(nfile, fid, beta_tag, 0, f_abMO_dic[0], float(FL_homo_ene[fid][0]))
        else:           # 0 for alpha, 1 for beta
            mplt_ab_draw(nfile, fid, beta_tag, 0, f_abMO_dic[0], float(FL_homo_ene[fid][0]))
            mplt_ab_draw(nfile, fid, beta_tag, 1, f_abMO_dic[1], float(FL_homo_ene[fid][1]))
            print f_abMO_dic[0], float(FL_homo_ene[fid][0])
            print f_abMO_dic[1], float(FL_homo_ene[fid][1])

        fid+=1

    #### draw MPL link line
    if L_link:
        print "Draw Link in ", os.path.basename(__file__)
        #### for each range in MO_link_hl_ene
        for n in range(nfile-1):
            nth_links=MO_link_hl_ene[n]
            x=Xrange_nf(nfile, n, 0, 0, 1)
            xl=x[1]     # 1=right side of 0-th file is left-link
            x=Xrange_nf(nfile, n+1, 0, 0, 1)
            xr=x[0]     # 0=left side of 1-th file is right-link 
            x_link=[xl, xr]
            for ene_pair in nth_links:
                mplt_line_link(x_link, ene_pair)

    plt.show()
    print "normal stop"	


def main():
    parser = argparse.ArgumentParser(description="to analyize QChem outfile and draw MO levels with several files")
    parser.add_argument("-f", "--outfiles", nargs="+", help="more than 1 out files")
    #parser.add_argument("-co", "--mocoeff", action='store_true', help="do coeff analysis")
    parser.add_argument("-ag", "--atomgroup", nargs='*', help="atom list of 2 groups for mo and link")
    parser.add_argument("-a", "--atomlists", nargs='*', help="atom lists with id for all files")
    parser.add_argument("-t", "--motypes", nargs='*', choices=['ALL','SEL','SUB','OR','AND','AND2'], help="types how to combine the selected atoms also in mp_file_anal.py")
    parser.add_argument("-bt", "--basetype", nargs='*', help="add atom basis")
    parser.add_argument("-l", "--link", action='store_true', help="make links for multiple files")
    parser.add_argument("-lt", "--linktypes", nargs='+', choices=['flow', 'lsplit', 'rsplit'], help="for the same molecule (flow) and combination with other molecule (split)")
    parser.add_argument("-sg", "--splitgroup", nargs="+", help="divide atom groups into two group for split link")
    args=parser.parse_args()

    print args.atomlists

    analyze_files(args.outfiles, args.atomlists, args. atomgroup, args.motypes , args.basetype, args.link, args.linktypes, args.splitgroup)
    return(0)



if __name__ == "__main__":
    main()
