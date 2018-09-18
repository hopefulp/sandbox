#!/home/jackjack5/epd/bin/python
### !/Users/JoonhoPark/anaconda/bin/python
### made by J. Park
### 2016. 10. 25 developed for radical
### 2016. 10. 28 developed for line
### 2016. 11. 14 developed for homo-reference
### 2016. 11. 15 developed for nbo BD
### 2016. 11. 21 developed for MO Coefficients analysis
### 2016. 11. 24 developing for nbo charge
### 2016. 11. 24 developing for MO linkage
### 2016. 12. 13 developing using Hash

import sys
import re
from mplt_qcdraw import *

#### read qchem output file upto 3
narg=len(sys.argv)
if narg <= 1 :
    print "Error with no arguments: input qchem outfile"
    print "Usage:: sys.argv[0] qchem.out (with job=sp)"
    exit(1)
FL_qchem_outf=[]
for file in sys.argv:
    if re.search("out", file):
    	FL_qchem_outf.append(file)
print "input files: ", FL_qchem_outf

######## analysis options
MOene=1
NBO_anal= 1             
MOCoeff_anal=1          
                        
######## MPL module mplt_qcdraw constants
####     FL stands for variable hase file list [ [f1:...],[f2:...],[f3:...] ]
FL_n_lumo=[1,1,2]        # cut number of lumo but include degeneracy of lumo of CO2
nfile=len(FL_qchem_outf)

'''
bonding with O
bonding                                             anti-bonding
1                                       88 homo    
2   92(w) 93(w) 94(s,bo,O1) 99(s,bb,C-O1)           100(s,bb-ab,C-O1
3   CO2                                             lumo

homo-1=-1 homo=0, lumo=1, lumo-1=2
MO_link_hl_id=((id_1st_mol, id_2nd_mol),(id_2nd_mol, id_3rd_mol))
'''
MO_link_id=(94,99)
MO_link_hl_id=[[[0,0],[0,1]], [[0,1],[1,1]]]
FL_sMO_dic=[]
FL_link_dic=[]
######## Key Word for Detection
####     key word for MO
#KW_MOene_occ_alpha="Alpha MOs"
KW_MOene_occ_alpha="Occupied"
KW_MOene_vir="Virtual"
KW_MOene_beta="Beta MOs"
#key_beta1="unrestricted"
KW_MOene_end="-----"

FL_homo_ene=[]      # it can have beta, save homo energy [ [f1 a, f1 b],[f2 a, f2 b]... ]
FL_homo_id=[]       # ene and MO id go parallel
FL_MOene_list=[]      # [ [e1, e2... ], [e1, e2...] [...
FL_MOene_list2d=[]
FL_MOene_select=[]
FL_MOid_select=[]
FL_beta_tag=[]      # 0 for restricted, 1 for beta spin of  unrestricted calculation 

######## KW for MO Coefficient
KW_MOcoeff="MOLECULAR ORBITAL COEFFICIENTS"
KW_nbasis="basis functions"
MOcoeff_crit=0.15

######## get id of MO link from NBO 
#### pair between 1-2 and 2-3 fragments
FL_Atom_types=[
    [['Ni','d']],           # 1st type
    [['C','p'],['O','p']]   # 2nd type
    ]
#### MO-kind is given by molecule section, just decide mol types
FL_MO_types=[]
for fname in FL_qchem_outf:
    if re.search("P", fname):   #### for P-P-N structure,   want to see Ni MO
        mo_type=0
    elif re.search("C", fname): #### for CO2                want to see C and O
        mo_type=1
    elif re.search("Ni", fname):
        mo_type=-1
    else:
        "Error for Molecular type for selective MO drawing"
        #exit(0)
    FL_MO_types.append(mo_type)
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

print_level=0
print_mo=1
"""
readline file and extract (1) energies, (2) coefficients, (3) nbos
"""
fid=0
# read three files
for file in FL_qchem_outf:
    tag_MO="YES"    # MO ene
    flag_MO="NO"
    flag_vline="NO"
    tag_beta="NO"
    tag_nbo="YES"
    flag_nbo="NO"
    flag_bd="NO"
    tag_moc="YES"
    flag_moc="NO"
    flag_moc_beta="OFF"
    inf=open(file, 'r')
    f_homo_i_line=[]
    f_list_ene_select=[]
    f_list_imo_select=[]
    f_ene_dic={}
    disp_MOids=[]
    line_ene_a=[]     # for job=optimization, whenever find new keyword, empty
    line_ene_b=[]
    i=0
    i_ene_line=0
    ibd=0
    ibd_atom=0
    nb_line=[]
    imoc=0
    imoc_block=0
    nbasis=0
    # read one file
    while 1:
        line=inf.readline()
        i += 1
        if not line: break         # end of while-loop, one file
        #### find nbasis for MO Coefficient here
        if MOCoeff_anal==1 and nbasis==0:
            if re.search(KW_nbasis, line) and re.search("shells", line):
                field=re.split("\s+", line)
                field=[x for x in field if x]
                nbasis=float(field[5])
                #print "found num of basis: ", nbasis
        ###### 1st: MO energies                         
        ################################ BLOCK 1 ###############################################
        #### find KW_MO energy block
        # inside one file
        if flag_MO == "NO" and tag_MO == "YES":
            if re.search(KW_MOene_occ_alpha, line):
                print "################### ", FL_qchem_outf[fid], "### BLOCK 1 :: QChem MO found ##################"
                flag_MO="YES"
                FL_beta_tag.append(0)
                eline_homo_tmp=[]
        #### have found MO
        elif flag_MO == "YES" and tag_MO == "YES":
            '''
            after found occupied orbital energy
            options except energies
            1. detect endline: change tag_MO into NO then this block is skipped
            2. virtual line
            '''
            #### if arrived the END line of MO energy, break
            if re.search(KW_MOene_end, line):
                tag_MO="Done"
            #### normally if found VIRTUAL, go to next line,
            elif re.search(KW_MOene_vir, line):
                flag_vline="YES"    # goes to the next line
            #### in case of Beta, there is blank line, how to treat BLANK
            #### have found beta
            elif not line.strip():
                flag_vline="NO"     # virtual energy ended 
            #### have found BETA, the same as above line
            elif re.search(KW_MOene_beta, line):
                tag_beta="YES"
                FL_beta_tag[fid]=1  # BETA is defined here 
                i_ene_line=0
                f_homo_i_line.append(iline_homo)
                print "spin beta in ", FL_qchem_outf[fid]
            #### save energy line
            else:
                # in case of beta, pass KW=="Occupied" again
                if re.search(KW_MOene_occ_alpha, line):
                    pass
                # save energy line before end line
                else:
                    if tag_beta == "YES":
                        line_ene_b.append(line)
                        if flag_vline != "YES":
                            iline_homo=i_ene_line
                    else:
                        line_ene_a.append(line)
                        if flag_vline != "YES":
                            iline_homo=i_ene_line
                            #print "alpha homo", iline_homo
                    #print i_ene_line, "here"
                    i_ene_line += 1
        ###### when MO is done, come to here
        ###### MO Coeff                         
        ################################ BLOCK 2 ########################################
        ###### select MOs
        if MOCoeff_anal==1 and tag_moc=="YES": # "cut" means alpha has done
            #### Look for KW_MO coeff
            if flag_moc == "NO" :
                if re.search(KW_MOcoeff, line) :
                    flag_moc="YES"
                    moc_block=[]
                    fl_homo_id=FL_homo_id[fid][0]
                    if flag_moc_beta == "ON":
                        fl_homo_id=FL_homo_id[fid][1]
                    print "################ BLOCK 2 :: MO Coeff starts "
                    #print FL_qchem_outf[fid], ":: MO Levels"
                    continue
            elif flag_moc == "YES":
                #### end of MO Coefficients
                if line =='\n':
                    if FL_beta_tag[fid]==0:
                        tag_moc="Done"
                    #elif flag_moc_beta=="OFF":
                    #    flag_moc="NO"
                    #    flag_moc_beta="ON"
                    #else:
                    #    tag_moc="Done"
                    continue
                #### save MO id
                if imoc==0:
                    #print line, 
                    imo=re.split("[\s\(\)-]+",line)
                    imo=[x for x in imo if x]
                    n_ene_col=len(imo)
                    imoc+=1
                    continue
                #### save energy
                elif imoc==1:
                    imo_ene=re.split("\s+",line)
                    imo_ene=[x for x in imo_ene if x]
                    del imo_ene[0]
                    #print imo_ene
                    imoc+=1
                    continue
                ###### arrive MOCoefficients block
                elif imoc <= 2+nbasis:
                    moc_line=re.split("[\s\(\)-]+", line) # coeff's are positive now
                    moc_line=[x for x in moc_line if x]
                    #print moc_line
                    # not storing but display
                    istart = len(moc_line) - n_ene_col 
                    for j in range(n_ene_col):
                        mo_symbol=[]
                        moc_tmp=[]
                        if float(moc_line[istart+j]) > MOcoeff_crit:
                            #print moc_line[1], moc_line[2],
                            mo_symbol.extend((moc_line[1], moc_line[2]))
                            if istart==3:
                                if print_level >= 1:
                                    print imo[j], imo_ene[j], moc_line[istart+j]
                                moc_tmp.extend((mo_symbol,imo[j],imo_ene[j],moc_line[istart+j]))
                            elif istart==4:
                                if print_level >= 1:
                                    print moc_line[3], imo[j], moc_line[istart+j]
                                mo_symbol.append(moc_line[3])
                                moc_tmp.extend((mo_symbol,imo[j],imo_ene[j],moc_line[istart+j]))
                            moc_block.append(moc_tmp)
                            #print moc_block
                    #### if read nbasis + 2 line, treat block and initialize
                    if imoc == 1+nbasis:
                        imoc=0
                        imoc_block+=1
                        #### treat block of MOC here
                        #### FORMAT of moc_block: [[P, 1, s], mo_id, mo_ene, mo_coeff
                        sorted_moc_block=sorted(moc_block, key=lambda x:float(x[1]))
                        for item in sorted_moc_block:
                            if print_level >=1:
                                print item
                        k=0     #### index for many coeff with the same MO id
                        s='-'
                        null="  "
                        id_tmp=0
                        
                        #### SELECT ENERGY LEVEL FOR MPL
                        for eff_coeff in sorted_moc_block:
                            #print "format:: [[atom_kind id angular_mom] mo_id energy MO_Coeff]"
                            #print eff_coeff
                            ######## Filter 1: by FL_n_lumo 
                            #### cut valence orbital
                            if int(eff_coeff[1]) > fl_homo_id+FL_n_lumo[fid]:
                                if FL_beta_tag[fid]==0:
                                    tag_moc="Done"
                                elif FL_beta_tag[fid]==1 and flag_moc_beta=="ON":
                                    tag_moc="Done"
                                #### run beta again
                                elif FL_beta_tag[fid]==1 and flag_moc_beta=="OFF":
                                    tag_moc="YES"
                                    flag_moc_beta="ON"
                                    flag_moc="NO"       # find KW for beta
                                break
                            ######## Filter 2: by atom-kind
                            #### cut if it's not Ni
                            if FL_MO_types[fid] == 0:
                                if eff_coeff[0][0] == 'Ni' or (eff_coeff[0][0]=='C' and eff_coeff[0][1]=="1") or (eff_coeff[0][0]=='O' and eff_coeff[0][1]=="2"):
                                    pass
                                else:
                                    continue
                                #elif eff_coeff[0][0]=='Ni' and eff_coeff[0][1]=='s':
                                #    continue
                            #### cut if it's not p-orbital for CO2
                            elif FL_MO_types[fid] == 1:
                                if len(eff_coeff[0])==2 and  eff_coeff[0][1]=='s':
                                    continue
                                if len(eff_coeff[0])==3 and eff_coeff[0][2]=='s':
                                    continue
                            else:
                                pass
                            #print eff_coeff
                            imo=eff_coeff[1]
                            mo_ene=eff_coeff[2]
                            #disp_MOids.append(id_mo)
                            if imo == id_tmp:
                                print "\t", s.join(eff_coeff[0]), "\t", eff_coeff[3]
                            else:
                                k=0
                                id_tmp=imo
                                f_list_imo_select.append(imo)
                                f_list_ene_select.append(mo_ene)
                                if print_mo==1:
                                    print "MO level: ", imo, mo_ene
                                    print "\t", s.join(eff_coeff[0]), "\t", eff_coeff[3]
                            #if tag_moc=="cut":
                            #    continue
                            k+=1
                        #print 
                        moc_block=[]
                        continue
                imoc+=1
                
        ###### Block 3 :: NBO BonDing extraction
        ######################### BLOCK 3 #####################################
        if NBO_anal==1:
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
        ############################ BLOCK 4 #################################################
        ###### Block 4 :: extract energy from MO energy lines
        ######              Homo Energy and ID
        if tag_MO=="Done":
            #print "######################### BLOCK 4 ###############################"
            #print FL_qchem_outf[fid], "MO obtain HOMO Energy and id"
            f_homo_i_line.append(iline_homo)
            #print f_homo_i_line
            if print_level >= 2:
                print FL_qchem_outf[fid], "homo energy line number: ", f_homo_i_line[fid]
            #### f is for 1 file, 2d means alpha, beta
            f_ene_line2d=[]
            f_ene_line2d.append(line_ene_a)
            if tag_beta == "YES":
                f_ene_line2d.append(line_ene_b)
            if print_level >= 2:
                print f_ene_line2d
            
            f_ene_tmp=[]

            list_ene=[]
            list_ene2d=[]   # 2d for [[alpha], [beta]]
            e_homo_tmp=[]
            e_homo_id_tmp=[]
            iab=0
            # split alpha beta
            for ab_line_list in f_ene_line2d:     # alpha beta [ A[line...],B[line...]]
                list_ene_tmp=[]
                j=0     # line index
                ihomo=0
                for e_line in ab_line_list:     # line index, it might be list[0] / [ list[0], list[1]]
                    #print j, e_line
                    ene=re.split("\s+", e_line)
                    one_line_ene=[]
                    for x in ene:
                        if x:                   # remove empty entry
                            list_ene.append(x)
                            one_line_ene.append(x)
                            list_ene_tmp.append(x)
                    #### finding homo ene level
                    if j == f_homo_i_line[iab]:
                        ehomo=list_ene[len(list_ene)-1]
                        ihomo=j*8+len(one_line_ene)
                        e_homo_tmp.append(ehomo)
                        e_homo_id_tmp.append(ihomo)
                        #print j, ehomo, ihomo
                    j+=1
                list_ene2d.append(list_ene_tmp)
                iab += 1
            FL_homo_ene.append(e_homo_tmp)
            FL_homo_id.append(e_homo_id_tmp)
            FL_MOene_list.append(list_ene)
            FL_MOene_list2d.append(list_ene2d)  # 2d is done in block 1 and 4
            print "HOMO energy, ID, beta tag:: ", FL_homo_ene[fid], FL_homo_id[fid], FL_beta_tag[fid]
            if 1 in  FL_beta_tag:
                print "Orbital energies in ", FL_qchem_outf[fid], ":: alpha, beta are combined"
            #### first block and last 4-th block are skipped
            tag_MO="Finish" 

        # end of last if-block
    # end of read 1 file
    inf.close() 
######## Filter 3: retailing plot list in each file
#### cut lower energy limit for drawing in figure
    del_id=[]
    k=0
    print f_list_imo_select, f_list_ene_select
    #exit(10)
    for imo, moe in zip(f_list_imo_select,f_list_ene_select) :
        if float(moe) < YMIN:
            del_id.append(k)
        k+=1
    #print del_id
    for id in reversed(del_id):
        del f_list_imo_select[id]
        del f_list_ene_select[id]
    #### make dictionary
    if FL_beta_tag[fid] == 0:
        for imo, moe in zip(f_list_imo_select,f_list_ene_select) :
            i_imo=int(imo)
            f_mo=float(moe)
            print i_imo, moe
            f_ene_dic[i_imo]=float(f_mo)
    else:
        #### change key for alpha, beta mo-id
        print "beta is not developed for dictionary"
        exit(20)
    print f_ene_dic
    FL_MOid_select.append(f_list_imo_select)
    FL_MOene_select.append(f_list_ene_select)
    FL_sMO_dic.append(f_ene_dic)
    fid += 1
######################################################
###### end of reading all  the files
print "############### Analized all the files ##########################"
print "Total:: HOMO energy, ID, beta tag: ", FL_homo_ene, FL_homo_id , FL_beta_tag

################# Block 5 ::PLOT using MPL ###########
fid=0
str="--------"
"""
MO_link_id=(94,99)
                    for left-side link      for right-side link
                 1st-mol 2nd-mol            2nd-mol 3rd-mol
MO_link_hl_id=[[[ 0       ,0  ],[0 ,   1]], [[0,       1],[    1,     1]]]
We need          y-val  y-val y-val y-val 
"""
print "Block5: Final MO energies for mpl plot:: "
from copy import copy, deepcopy
MO_link_hl_ene=deepcopy(MO_link_hl_id)

for list_ene, list_imo in zip(FL_MOene_select, FL_MOid_select) :
    print FL_qchem_outf[fid], "MO id and energy:"
    #### scan all the needed energy
    if nfile == 3:
        if fid==0:
            for k in range(len(MO_link_hl_id[0])):
                ####       0 for link side, k for several MO levels, 0 for 1st-mol, 
                print MO_link_hl_id[0][k][0]
        
    for m_id, m_ene in zip(list_imo, list_ene):
        print "%5d %10.3f\t" % (int(m_id), float(m_ene))
    #print 
    #print FL_qchem_outf[fid],"energy", list_ene
    #print FL_qchem_outf[fid], "mo id\t", str.join(list_imo)
    ### initialize variables for each file drawing
    tag_draw_beta="NO"
    xlength=0.4
    x_3kinds=[]
    x_3kinds_b=[]
    ### for single level
    if FL_beta_tag[fid] == 1:
        xlength /= 2
    x_l = XMIN + 0.2	+ X_SHIFT[fid]
    x_r = x_l  + xlength
    x_d = x_r - x_l	
    x_3kinds.append([x_l, x_r])
    if FL_beta_tag[fid] == 1:
	    x_3kinds_b.append([x_l+AB_GAP, x_l+AB_GAP+x_d])
    ### for double degeneracy
    x_l1 = x_l		    	
    x_r1 = x_l + xlength/3  
    x_d  = x_r1-x_l1
    x_l2 = x_r1 + xlength/3 
    x_r2 = x_l2+x_d		  
    x2l=[x_l1, x_r1]
    x2r=[x_l2, x_r2]
    x2=x2l+x2r
    x_3kinds.append(x2)
    #print "print x2 interval: ", x2
    if print_level=="ON":
        print "x2 range:: ", "%.2f "*len(x2) % tuple(x2)
    if FL_beta_tag[fid] == 1:
        xbeta2l= [x_l1+AB_GAP, x_r1+AB_GAP]
        xbeta2r= [x_l2+AB_GAP, x_r2+AB_GAP]
        x2beta=xbeta2l+xbeta2r
        x_3kinds_b.append(x2beta)
        if print_level=="ON":
            print "x2b range:: ", "%.2f "*len(x2beta) % tuple(x2beta)
    ### for triple degeneracy
    x_r1 = x_l + xlength/5
    x_d  = x_r1-x_l
    x_l2 = x_r1+0.05
    x_r2 = x_l2+x_d
    x_l3 = x_r2+0.05 
    x_r3 = x_l3+x_d
    x3l=[x_l1, x_r1]
    x3m=[x_l2, x_r2]
    x3r=[x_l3, x_r3]
    x3=x3l+x3m+x3r
    x_3kinds.append(x3)
    if print_level=="ON":
        print "x3 range:: ", "%.2f "*len(x3) % tuple(x3)
    if FL_beta_tag[fid] == 1:
        xbeta3l=[x_l1+AB_GAP, x_r1+AB_GAP]
        xbeta3m=[x_l2+AB_GAP, x_r2+AB_GAP]
        xbeta3r=[x_l3+AB_GAP, x_r3+AB_GAP]
        x3beta=xbeta3l+xbeta3m+xbeta3r
        x_3kinds_b.append(x3beta)
        if print_level=="ON":
            print "x3b range:: ", "%.2f "*len(x3beta) % tuple(x3beta)
    i=0
    multi_degen=1
    e_homo=float(FL_homo_ene[fid][0])
    ### for H radical
    if len(list_ene) == 1:
	    y=[float(list_ene[0]), float(list_ene[0])]
	    f_draw(1, x_3kinds[0], y, e_homo)
	    fid+=1
	    continue
    temp_ene=-1000
    tag_last_draw="NO"
    for v_ene in list_ene:
        if re.search("\*", v_ene): # in qchem outfile, 1st low number might be ******
            i+=1
            continue
        v_ene = float(v_ene)
    	if temp_ene==-1000:         # for initial energy value, just save 
            temp_ene=v_ene
            i+=1
            continue
    	y=[temp_ene, temp_ene]
    	##### if energies are same, increase multi degeneracy 
        if temp_ene == v_ene:
            multi_degen += 1
            ### if not last energy
            if i != len(list_ene)-1:
                i+=1
                continue
	        ### draw 3 plots if last ene and finish
            else:
                f_draw(multi_degen, x_3kinds[multi_degen-1], y, e_homo)
		    #print "finish drawing"
	    ### if energies are different 
        else:
            ### draw 3 kind of plots with degeneracy for previous energies
            if multi_degen > 3:
                print "so many degeneracy in drawing"
                exit(10)
            f_draw(multi_degen, x_3kinds[multi_degen-1], y, e_homo)
            # if last ene and no degenracy, plot additional energy
            if i == len(list_ene)-1:
                y=[v_ene, v_ene]
                f_draw(1, x_3kinds[0], y, e_homo )
            # change into beta draw
            if(temp_ene > v_ene):
                print "Draw beta spin"
                print temp_ene, v_ene
                tag_draw_beta="YES"
                e_homo=float(FL_homo_ene[fid][1])
                x_3kinds=x_3kinds_b
        i+=1	
    	temp_ene=v_ene
    	multi_degen = 1
        ### end of one energy value
    fid+=1
    ### end of energy list
    
plt.show()
print "normal stop"	

