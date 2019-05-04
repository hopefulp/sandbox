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
### required: draw only a few levels    
### saved in Tmppy/draw_MOs_link_dev?.py
### to control the display of MO levels, CTRL1, CTRL2
### git test
### argparse

import re
import os
import argparse
import _math
from common import *
from mplt_qcdraw import *
from mp_mo_dump import *
from mp_file_anal import *
from mplt_mo_ini import *
from atom_valence import *

#### print option
V_print_link=0
V_print = 0


def f_extract_ene(index_homo_line, line_ene_ab):
    """
        make list of MO_ene
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
        j=0     # line index
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
            #### finding homo ene level
            if j == index_homo_line[iab]:
                ehomo=list_ene[len(list_ene)-1]
                ihomo=j*Nene_line_QC+len(one_line_ene)
                e_homo_ab.append(ehomo)
                e_homo_id_ab.append(ihomo)
            j+=1
        # after finish alpha/beta            
        list_ene_ab.append(list_ene)     # store in 2D
        iab += 1
    # after finish a/b list        
    print "HOMO energy, ID, beta, max_lumo in function:", whoami(), e_homo_ab, e_homo_id_ab, beta
    return list_ene_ab, e_homo_ab, e_homo_id_ab


################### MAIN LOOP ####################
def analyze_files(files,atomlists, motypes ,L_link): 
    # start analyze_files function 
    print files, "MO filter[atom|type]", atomlists, motypes, "L_link = ", L_link 
   
    """
        FL_:: all-files variable [[[file1-a],[file1-b]],[[file2-a],[file2-b]],...
              w. a,b dependency,  last dimension is 2
        Fl_:: all-files variable [[file1],[file2],...
        fL_:: for single file variable with a, b dependency
        f_:: for single file variable 
        nlumo for plot is i_homo_a + Nmax_lumo
    """
    Fl_qchem_outf = files[:]
    nfile=len(Fl_qchem_outf)
    ### FILTER MO for selected atoms
    #   if atomlists are defined, MOC analysis is auto-activated
    if atomlists:                   # go to MOC BLOCK 2
        L_mo_coeff = True
        MOcoeff_crit=0.15           # 0.15, 0 for all the levels
    else:
        L_mo_coeff = False
        max_lumo=Nmax_virtual
    
    if L_link:
        MOlink_hl_id=MO_link_hl_id[:]

    FL_sMO_dic=[]
    FL_link_dic=[]

    ### MO block
    FL_norbital=[]
    FL_homo_ene=[]      # it can have beta, save homo energy [ [f1 a, f1 b],[f2 a, f2 b]... ]
    FL_homo_id=[]       # ene and MO id go parallel
    FL_MOene_list=[]      # [ [e1, e2... ], [e1, e2...] [...
    FL_MOene_select=[]
    FL_abMO_dic_list=[]
    FL_MOid_select=[]
    Fl_beta_tag=[]      # 0 for restricted, 1 for beta spin of  unrestricted calculation 

    ### MOC block
    Fl_MO_types=[]      # [ONE(f1), 1sub(f2)...
    Fl_MO_atoms=[]      # [f1: [Ni, P], f2: [C, O], ..
    Fl_nlumo=[]         # [nlumo-f1, nlumo-f2...    same as a|b

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

    """
    Constants in module
    readline file and extract (1) energies, (2) coefficients, (3) nbos
    """
    fid=0
    # read three files
    for file in Fl_qchem_outf:
        # tags to get MO energies
        tag_MO="YES"    # YES, DONE, FINISH MO ene
        flag_MO="OFF"
        # 
        tag_nbo="YES"
        flag_nbo="NO"
        flag_bd="NO"
        # tags to get MO coefficients
        tag_moc="YES"
        flag_moc="OFF"
        flag_moc_beta="OFF"
        f_homo_i_line=[]        # one or two(beta_tag) values of index of homo line
        f_list_ene_select=[]
        f_list_ene_sel_ab=[]
        f_list_imo_select=[]
        f_list_imo_sel_ab=[]
        f_ene_dic={}
        f_ene_dic_beta={}
        f_dic_ab_list=[]
        f_norbital_ab=[]
        disp_MOids=[]
        f_line_ene_a=[]     # for job=optimization, whenever find new keyword, empty
        f_line_ene_b=[]
        fL_line_ene=[]
        fL_homo_id=[]
        ibd=0
        ibd_atom=0
        nb_line=[]
        nbasis=0
        temp_id_atom='Z'    # not atom id such as Ni, Ni-1, Ni-2...
        #### for MOCoefficient of BLOCK 2
        imoc=0
        imoc_block=0
        id_tmp=0    # mo-id starts from 1 in MO Coefficient Block
        mo_ene_tmp=-100
        dump_mo_id_block=[]

        s='-'       # atom-id linker
        # read one file
        inf=open(file, 'r')
        i=0
        b=0     # beta ene line index
        while 1:
            line=inf.readline()
            i += 1
            if not line: break         # end of while-loop, one file
            #### find nbasis for MO Coefficient here
            if L_mo_coeff==1 and nbasis==0:
                if re.search(KW_nbasis, line) and re.search("shells", line):
                    field=re.split("\s+", line)
                    field=[x for x in field if x]
                    nbasis=float(field[5])
                    #print "found num of basis: ", nbasis
# func
    # for each file
        # while
            ###### 1st: MO energies                         
            ###### BLOCK 1 :: MO energy block
            # inside one file
            if flag_MO == "OFF" and tag_MO == "YES":
                ### found MO block (first "occupied" is for Alpha spin)
                if re.search(KW_MOene_occ, line):
                    print "################### ", Fl_qchem_outf[fid], "#########################################"
                    print "###################  BLOCK 1 :: QChem MO energy"
                    flag_MO="ON"
                    Fl_beta_tag.append(0)   # beta option is initiated
                    ### initiate all the variables used in this block
                    i_ene_line=0            # energy line index starts from 0
                    tag_beta="NO"
                    flag_vline="NO"         # now occupied orbital
                    continue
            #### have found MO
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
                    f_homo_i_line.append(iline_homo) # twiced if beta_tag, iline_homo is assigend below
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
                        if b == 0:
                            print "check beta line", line.rstrip, f_line_ene_a[b]
                            if line.rstrip() == f_line_ene_a[b]:
                                print "Beta spin is detected but same as Alpha spin"
                                Fl_beta_tag[fid]=0
                                tag_MO="Done"
                                continue
                        f_line_ene_b.append(line.rstrip())
                        if flag_vline != "YES":
                            iline_homo=i_ene_line
                        b+=1
                    else:
                        f_line_ene_a.append(line.rstrip())
                        if flag_vline != "YES":
                            iline_homo=i_ene_line
                    i_ene_line += 1
                    continue
            elif tag_MO == 'Done':
                ### in case that alpha and beta are same though restricted
                #if f_line_ene_a == f_line_ene_b:
                #    Fl_beta_tag[fid] = 0

                tag_MO = 'Finish'
                fL_line_ene.append(f_line_ene_a)
                if Fl_beta_tag[fid] == 1:
                    fL_line_ene.append(f_line_ene_b)
                fL_ene_ab,f_homo_ab,f_homo_id_ab=f_extract_ene(f_homo_i_line, fL_line_ene)
                FL_homo_ene.append(f_homo_ab)
                FL_homo_id.append(f_homo_id_ab)
                FL_MOene_list.append(fL_ene_ab)
                
                f_norbital_ab.append(len(fL_ene_ab[0]))
                if Fl_beta_tag[fid] == 1:
                    f_norbital_ab.append(len(fL_ene_ab[1]))
                FL_norbital.append(f_norbital_ab)
                continue
            ################################ BLOCK 2 :: MO Coeff ########################################
            ###### ene & index are obtained here
            if L_mo_coeff==1: # "cut" means alpha has done
                #### Look for KW_MO coeff
                if flag_moc == "OFF" and  tag_moc=="YES":
                    if re.search(KW_MOcoeff, line) :    #common for alpha & beta
                        ### Initiation for MOC
                        max_nlumo = Nmax_virtual
                        Fl_nlumo.append(max_nlumo)
                        #### obtain atoms for MO and mol types in QChem outfile
                        if atomlists:
                            fl_mo_atoms=atomlists[fid].split('\s+')
                            Fl_MO_atoms.append(fl_mo_atoms)
                        if motypes:
                            f_mo_type=motypes[fid]
                        else:
                            f_mo_type=get_atom2motype(fl_mo_atoms)         #Fl_MO_atoms[:]
                        Fl_MO_types.append(f_mo_type)             #Fl_MO_type[:]
                        #Fl_MO_atoms, Fl_MO_type = get_mo_labels(files, atoms, i_mo_types)
                        Fl_nlumo, MO_link_hl_id = get_nlumo_linkid(nfile)

                        flag_moc="ON"
                        moc_block=[]
                        ### use common homo index for a|b
                        fL_homo_id=FL_homo_id[fid]      # list of one|two elements
                        f_homo_max = FL_homo_id[fid][0] # Alpha has always large value
                        f_lumo_id_cut=f_homo_max+max_nlumo
                        #lumo_id_cut=f_homo_id+Fl_nlumo[fid]
                        if f_norbital_ab <= f_lumo_id_cut:
                            print "Error:: too small lumos"
                            exit(33)
                        print "MO cut id: ", fL_homo_id, f_lumo_id_cut
                        if flag_moc_beta == "ON":
                            f_homo_id=FL_homo_id[fid][1]
                            print "###################  BLOCK 2 :: MO Coeff beta starts "
                            imoc=0
                            imoc_block=0
                            id_tmp=0            # mo-id starts from 1 in MO Coefficient Block
                            mo_ene_tmp=-100
                            dump_mo_id_block=[]
                            f_list_imo_sel_ab=[]
                            f_list_ene_sel_ab=[]
                        else:
                            print "###################  BLOCK 2 :: MO Coeff starts "
                            print "MO atoms|type", fl_mo_atoms,"|",f_mo_type
                        continue
                elif flag_moc == "ON":
                    #### if end of MO Coefficients
                    if line =='\n':
                        if Fl_beta_tag[fid]==0:
                            tag_moc="Done"
                        continue
                    #### 1. new MOC block and 1st line: 1.1 save MO id
                    if imoc==0:
                        #imo=re.split("[\s\(\)]+",line)     # for parenthesis
                        imo=re.split("[\s]+",line)
                        imo=[x for x in imo if x]
                        n_ene_col=len(imo)
                        #print "imoc==0; n_ele_col=", n_ene_col, line 
                    #### save energy
                    elif imoc==1:
                        imo_ene=re.split("\s+",line)
                        imo_ene=[x for x in imo_ene if x]
                        del imo_ene[0]
                        #print imo_ene
                    ###### arrived MO Coefficient line in MOC block: read moc line 
                    elif imoc <= 2+nbasis:
                        #moc_line=re.split("[\s\(\)-]+", line) # coeff's are positive now
                        moc_line=re.split("[\s]+", line)    # coeff's might be negative 
                        moc_line=[x for x in moc_line if x]
                        #print moc_line
                        #not storing but display
                        istart = len(moc_line) - n_ene_col 
                        id_atom=moc_line[1]
                        if re.match("\d",moc_line[2]):
                            id_atom+=moc_line[2]        # id_atom == Ni1, Ni2, H30, H31, etc
                        #### FILTER 1  skip core orbital based on atom species in moc_line[1]: follow QChem basis
                        if id_atom != temp_id_atom:
                            temp_id_atom=id_atom
                            ic_atom=0
                            print moc_line
                            ncore=Atom_Core_631g[Atom_Table.index(moc_line[1])]
                            print ncore
                            #### if there is beta electron, makes error
                        else:
                            ic_atom+=1
                        #print "atom ID & num of core atomic orbital: ", id_atom, ncore
                        if ic_atom <= ncore:
                            pass
                        else:
                            #### for each MOC line, SAVE n_ene_col==6 (if last, n_ene_col<=6) in moc_block
                            for j in range(n_ene_col):
                                mo_symbol=[]
                                moc_tmp=[]
                                #### FILTER 2 of MOC line by MOcoeff_crit
                                if abs(float(moc_line[istart+j])) > MOcoeff_crit:
                                    mo_symbol.extend((moc_line[1], moc_line[2]))
                                    if istart==3:
                                        #if V_print >= 3:
                                        #    print imo[j], imo_ene[j], moc_line[istart+j]
                                        moc_tmp.extend((mo_symbol,imo[j],imo_ene[j],moc_line[istart+j]))
                                    elif istart==4:
                                        #if V_print >= 3:
                                        #    print moc_line[3], imo[j], moc_line[istart+j]
                                        mo_symbol.append(moc_line[3])
                                        moc_tmp.extend((mo_symbol,imo[j],imo_ene[j],moc_line[istart+j]))
                                    elif istart==5:
                                        #if V_print >= 3:
                                        #    print moc_line[3], imo[j], moc_line[istart+j]
                                        mo_symbol.extend([moc_line[3],moc_line[4]])
                                        moc_tmp.extend((mo_symbol,imo[j],imo_ene[j],moc_line[istart+j]))
                                    else:
                                        print "MO Coefficients format error"
                                        exit(11)
                                    #print mo_symbol,
                                    #### Add moc line in moc_block of 6 mo_id
                                    moc_block.append(moc_tmp)
                                    #print "1st block in MO Coeff"
                                    #print(np.matrix(moc_block))
# func
    # for file
        # while
            #if MOC true
                # if moc_flag
                    # if
                        #### if read nbasis + 2 line
                        #### HAVE READ one MOC block with 6 mo_id and TREAT it and INITIALIZE
                        if imoc == 1+nbasis:
                            imoc=0
                            imoc_block+=1
                            #### treat block of MOC here
                            #### FORMAT of moc_block: [ [   P,     1,           s], mo_id, mo_ene, mo_coeff
                            ####                      [ [atom, index, angular mom], mo_id, mo_ene, mo_coeff
                            ####                                         SORTED by  mo_id
                            #print moc_block
                            sorted_moc_block=sorted(moc_block, key=lambda x:float(x[1]))
                            for item in sorted_moc_block:
                                if V_print >=2:
                                    print item # [['Ni', 'd', '4'], '100', '-0.016', '0.19560']
                            
                            #### eff_coeff for effective coefficients == filtered by valuable coefficient
                            #### SORTED_MOC_BLOCK with valence shell, effective coefficients
                            for eff_coeff in sorted_moc_block:
                                #print "format:: [[atom_kind id angular_mom] mo_id energy MO_Coeff]"
                                #print eff_coeff
                                ######## FILTER 3: cut lines using FL_nlumo[0] 
                                #print int(eff_coeff[1]) ,"<", lumo_id_cut
                                if int(eff_coeff[1]) > f_lumo_id_cut:
                                    #### dump MO Coefficient before break
                                    #print "last: ", f_list_imo_sel_ab
                                    imo_select=f_imo_dump(dump_mo_id_block, Fl_MO_atoms[fid], Fl_MO_types[fid],id_tmp)
                                    if imo_select:
                                        f_list_imo_sel_ab.append(id_tmp)
                                        f_list_ene_sel_ab.append(mo_ene_tmp)
                                        #print "last-1: ", f_list_imo_sel_ab
                                    #### save imo & ene for alpha or beta
                                    f_list_imo_select.append(f_list_imo_sel_ab)
                                    f_list_ene_select.append(f_list_ene_sel_ab)
                                    #print f_list_imo_select[0]
                                    if Fl_beta_tag[fid]==0:
                                        tag_moc="Done"
                                    elif Fl_beta_tag[fid]==1 and flag_moc_beta=="ON":
                                        tag_moc="Done"
                                    #### run beta again
                                    elif Fl_beta_tag[fid]==1 and flag_moc_beta=="OFF":
                                        tag_moc="YES"
                                        flag_moc_beta="ON"
                                        flag_moc="OFF"       # find KW for beta
                                    #### FINISH MO Coeff selection
                                    break
                                ######## DUMP for selected atom line in selected 6 mo-id's ############################
                                ######## FILTER 4       :: cut by atom species here
                                ####     f_imo_dump()   :: decision for dumping the imo-block
                                #### cf. eff_coeff== [['Ni', 'd', '4'], '100', '-0.016', '0.19560']
                                #### 'pass' writes the line in dump string      of str_dump
                                #### 'continue' skips the line before writing   in str_dump
                                if Fl_MO_types[fid] == "ONE":
                                    if eff_coeff[0][0] == Fl_MO_atoms[fid][0]:
                                        pass
                                    else:
                                        continue
                                #### if eff_coeff[0][0] is one of atom_list, save
                                #### 'ALL' and '1sub' works the same but different in f_imo_dump()
                                #### more detail screening difference between 'ALL' and '1sub', refer to f_imo_dump
                                elif Fl_MO_types[fid] == "ALL" or Fl_MO_types[fid] == "1sub":    # 'all' for the selected atom 
                                    f_tag=0
                                    for atom in Fl_MO_atoms[fid]:
                                        if eff_coeff[0][0] == atom:
                                            f_tag=1
                                            break
                                    if f_tag==1:
                                        pass
                                    else:
                                        continue
                                #### for more options for atom id w. index
                                elif Fl_MO_types[fid] == "1sub-a":
                                    if (eff_coeff[0][0] == Fl_MO_atoms[fid][0]) or (eff_coeff[0][0]==Fl_MO_atoms[fid][1] and eff_coeff[0][1]=="1") or (eff_coeff[0][0]==Fl_MO_atoms[fid][2]):
                                        pass
                                    else:
                                        #if eff_coeff[0][0]=='Ni' and eff_coeff[0][1]=='s':
                                        continue    # inside inner for-eff_coeff
                                elif Fl_MO_types[fid] == 'NONE': # for CO2
                                    pass
                                else:
                                    "Error:: there is no atom type"
                                    exit(11)
                                #print eff_coeff
                                imo=int(eff_coeff[1])
                                mo_ene=eff_coeff[2]

                                #### if index MO is the same, save for write
                                if imo == id_tmp:
                                    if float(eff_coeff[3]) < 0:
                                        dump_str="      "+s.join(eff_coeff[0])+"\t    "+eff_coeff[3]
                                    else: 
                                        dump_str="      "+s.join(eff_coeff[0])+"\t     "+eff_coeff[3]
                                    dump_mo_id_block.append(dump_str)
                                #### if index MO is changed, dump
                                else:
                                    #### dump the previous saved MO Coefficients
                                    #### FILTER 4.5 :: the imo-block will be dumped?
                                    #### imo_select returns 1 if dump or 0 for skip
                                    imo_select=f_imo_dump(dump_mo_id_block, Fl_MO_atoms[fid], Fl_MO_types[fid], id_tmp)
                                    #### RECORD and INITIALIZE
                                    #print  id_tmp, mo_ene_tmp, imo_sel
                                    #### imo_select was decided to be printed so save
                                    if imo_select:
                                        #print id_tmp
                                        f_list_imo_sel_ab.append(id_tmp)
                                        f_list_ene_sel_ab.append(mo_ene_tmp)
                                        #print f_list_imo_sel_ab
                                    #### initialize mo id-1-block
                                    dump_mo_id_block=[]
                                    id_tmp=imo
                                    mo_ene_tmp=float(mo_ene)
                                    #### selection position 1
                                    #f_list_imo_sel_ab.append(imo)
                                    #f_list_ene_sel_ab.append(mo_ene)
                                    print imo       # str(imo) makes error ?
                                    dump_str="MO level  "+ "   "+mo_ene
                                    dump_mo_id_block.append(dump_str)
                                    if float(eff_coeff[3]) < 0:
                                        dump_str="      "+s.join(eff_coeff[0])+"\t    "+eff_coeff[3]
                                    else:
                                        dump_str="      "+s.join(eff_coeff[0])+"\t     "+eff_coeff[3]
                                    dump_mo_id_block.append(dump_str)

                            #### end of MO coeff one block
                            moc_block=[]
                            continue
                        #### have read one block of nbasis lines
                    #### inside MOC block
                    imoc+=1
                    continue
# func
    # for file
        # while
            # if
            ###### Block 3 :: NBO BonDing extraction
            ######################### BLOCK 3 #####################################
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
        #### change into all levels for T_select_bond=='all'
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
            
                f_list_ene_select.append(tmp_elist[:])
                f_list_imo_select.append(tmp_list[:])
        if V_print >= 3:
            print "After imo selection, make dictionary: ", f_list_imo_select
            print f_list_ene_select
        print "\"Draw MO levels of ", Fl_qchem_outf[fid], "\""
        print "MO_level_ID   Energy"
        #print T_select_bond
        for imo, moe in zip(f_list_imo_select[0],f_list_ene_select[0]) :
            i_imo=imo
            f_mo=moe
            print "%11d %8.3f" % (i_imo, moe)
            f_ene_dic[i_imo]=f_mo
        f_dic_ab_list.append(f_ene_dic)
        if Fl_beta_tag[fid] == 1:
            #### change key for alpha, beta mo-id
            print f_list_imo_select[1]
            print f_list_ene_select[1], len(f_list_ene_select[1])
            for imo, moe in zip(f_list_imo_select[1],f_list_ene_select[1]) :
                i_imo=imo
                f_mo=moe
                print i_imo, moe
                f_ene_dic_beta[i_imo]=f_mo
            f_dic_ab_list.append(f_ene_dic_beta)
        print f_dic_ab_list            
        #### code for beta is developed up to here
        #print f_ene_dic
        FL_MOid_select.append(f_list_imo_select[0])
        FL_MOene_select.append(f_list_ene_select[0])
        FL_sMO_dic.append(f_ene_dic)
        FL_abMO_dic_list.append(f_dic_ab_list)
        fid += 1
    ######################################################
    ###### end of reading all  the files
    print "############################################"
    print "############### Result #####################"
    print "        Filename        HOMO_energy & ID  beta_tag "
    for i in range(nfile): 
        print "%20s" % Fl_qchem_outf[i], 
        print FL_homo_ene[i], FL_homo_id[i], Fl_beta_tag[i]

    ################# Block 5 ::PLOT using MPL ###########
    str="--------"
    """
    MO_link_id=(94,99)
                        for left-side link      for right-side link
                     1st-mol 2nd-mol            2nd-mol 3rd-mol
    MOlink_hl_id=[[[ 0       ,0  ],[0 ,   1]], [[0,       1],[    1,     1]]]
    We need          y-val  y-val y-val y-val 
    """
    print "##### Block5: Final MO energies for mpl plot:: "


    """
    Link is only working for non-degenerate MOs
    """
    MO_link_hl_ene=[]

    ######## Make MO-link pair
    ######## how to obtain link pair

    #### To use mpl-plot, it needs energies for y-axis, not id
    if L_link and nfile >=2:
        if 'MOlink_hl_id' not in locals():
            print "Error:: No link pair input in", os.path_basename(__file__)
            exit(90)

        for i in range(nfile-1):        
            #iside_link=0     # 0 for between 0th and 1st files
            ihomo_f0=FL_homo_id[i][0]       # f0 for left file, f1 for right,
            ihomo_f1=FL_homo_id[i+1][0]     # due to alpha, beta
            #### here obtain link pair
            print "Link ID in %d-th range:: " % i  
            #### find MO id using homo-id
            #### between MO 0 and MO 2
            #### obtain link pair here
            #f_link
            link_y=[]       # make link pair
            for l_pair in MOlink_hl_id[i]:
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
        if V_print_L >= 0:
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
    parser.add_argument("-l", "--link", action='store_true', help="make links for multiple files")
    #parser.add_argument("-b", "--bonds", default='a', choices=['a','s'], help="selective")
    parser.add_argument("-co", "--mocoeff", action='store_true', help="do coeff analysis")
    #groupa = parser.add_mutually_exclusive_group()
    parser.add_argument("-a", "--atomlists", nargs='*', help="2D atom list for selective MO drawing")
    parser.add_argument("-t", "--motypes", nargs='*', help="types are written in mp_file_anal.py")
    #parser.add_argument("-t", "--motype", default='ONE', choices=['ONE','1sub','1sub-a','ALL','NONE'], help="MO-type depending on the atoms")
    #parser.add_argument("-t", "--motypes", default=0, type=int, nargs="+", help="number of motype should match with nfiles")
    args=parser.parse_args()
    '''
    if args.motypes:
        if not len(args.outfiles) == len(args.motypes):
            print "number of type should be equal to that of files\n"
            exit(1)
        FL_selection=args.motypes
    elif args.atomlists:
        if not len(args.outfiles) == len(args.atomlists):
            print "number of type should be equal to that of files\n"
            exit(1)
        FL_selection=args.atomlists
    '''

    #analyze_files(args.outfiles, args.mocoeff, args.atomlists, args.motypes , args.link)
    analyze_files(args.outfiles, args.atomlists, args.motypes , args.link)
    return(0)



if __name__ == "__main__":
    main()
