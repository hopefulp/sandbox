#!/home/jackjack5/epd/bin/python
### !/Users/JoonhoPark/anaconda/bin/python
### made by J. Park
### 2016. 10. 25 developed for radical
### 2016. 10. 28 developed for line
### 2016. 11. 14 developed for homo-reference
### 2016. 11. 15 developed for nbo BD
### 2016. 11. 21 developing for MO Coefficients analysis

import sys
import re
from mplt_qcdraw import *

#### read qchem output file upto 3
nfile=len(sys.argv)
#print sys.argv
if nfile <= 1 :
    print "Error with no arguments: input qchem outfile"
    print "Usage:: sys.argv[0] qchem.out (with job=sp)"
    exit(1)
list_qchem_ofiles=[]
for file in sys.argv:
    if re.search("out", file):
    	list_qchem_ofiles.append(file)
print "input files: ", list_qchem_ofiles

NBO_anal=1
MOCoeff_anal=1

#### MPL module mplt_qcdraw constants
nfile=len(list_qchem_ofiles)
#ini_xlimit(nfile)

#### global variable 
#keyword_a_occ="Alpha MOs"
keyword_a_occ="Occupied"
keyword_a_virt="Virtual"
keyword_beta="Beta MOs"
#key_beta1="unrestricted"
keyword_end="-----"
FL_homo_i_line=[]       # it can have beta
FL_homo_ene=[]          # it can have beta, save homo energy [ [f1 a, f1 b],[f2 a, f2 b]... ]
FL_ene_line=[]          # [ [eline1, eline2,...  ], [eline1, eline2...] [...
FL_list_ene=[]      # [ [e1, e2... ], [e1, e2...] [...
FL_ene_select=[]
FL_disp_MO_id=[]        
FL_disp_MO_ene=[]       
FL_beta_tag=[0,0,0]  # [1, 0, 1] 1 for beta spin of  unrestricted calculation 
### define orbital from MO coefficients
keyword_moc="MOLECULAR ORBITAL COEFFICIENTS"
keyword_nbasis="basis functions"
moc_crit=0.2

###### get id of MO link from NBO 
### pair between 1-2 and 2-3 fragments
FL_MO_types=[
    [['Ni','d']],           # 1st mol
    [['Ni','d']],           # 2nd mol
    [['C','p'],['O','p']]   # 3rd mol
    ]

nbo_atoms=['C', 'O']    # atom for 1st fragment, 2nd fragment

keyword_nbo="Bond orbital"
keyword_nbo_BD="BD"
keyword_nbo_end="NBO state"

AngSym=('1s', '2s', '2px', '2py', '2pz')

lcmo_mol=1      
lcmo_nlnh=[1,1]
lcmoA=[]
lcmoB=[]
lcmoC=[]

print_level=0
"""
readline file and extract (1) energies, (2) coefficients, (3) nbos
"""
fid=0
for file in list_qchem_ofiles:
    tag_MO="YES"
    flag_MO="NO"
    flag_vline="NO"
    tag_beta="NO"
    tag_nbo="YES"
    flag_nbo="NO"
    flag_bd="NO"
    tag_moc="YES"
    flag_moc="NO"
    inf=open(file, 'r')
    iline_homo_tmpl=[]
    f_list_ene_select=[]
    disp_MOids=[]
    line_ene=[]     # for job=optimization, whenever find new keyword, empty
    line_ene_b=[]
    i=0
    i_ene_line=0
    ibd=0
    imoc=0
    imoc_block=0
    nbasis=0
    moc_block=[]
    while 1:
        line=inf.readline()
        if not line: break         # end of while-loop, one file
        ### find nbasis for MO Coefficient here
        if MOCoeff_anal==1:
            if nbasis==0:
                if re.search(keyword_nbasis, line) and re.search("shells", line):
                    field=re.split("\s+", line)
                    field=[x for x in field if x]
                    nbasis=float(field[5])
                    #print "found num of basis: ", nbasis
        ###### 1st: MO energies
        ####################################################################################
        if flag_MO == "NO" and tag_MO == "YES":
            if re.search(keyword_a_occ, line):
                flag_MO="YES"
                eline_homo_tmp=[]
        ### found MO
        elif flag_MO == "YES" and tag_MO == "YES":
            # normally if found VIRTUAL, go to next line,
            if re.search(keyword_a_virt, line):
                flag_vline="YES"
            # how to treat BLANK
            elif not line.strip():
                pass
            # in case BETA
            elif re.search(keyword_beta, line):
                tag_beta="YES" 
                flag_vline="NO"   # for beta lumo
                FL_beta_tag[fid]=1
                i_ene_line=0
                iline_homo_tmpl.append(iline_homo)
                print "spin beta in ", list_qchem_ofiles[fid]
            # if arrived the END line of MO energy, break
            elif re.search(keyword_end, line):
                tag_MO="Done"
                #flag_nbo="YES"
                #break
            else:
                # in case of beta, pass Occupied again
                if re.search(keyword_a_occ, line):
                    pass
                # save energy line before end line
                else:
                    if tag_beta == "YES":
                        line_ene_b.append(line)
                        if flag_vline != "YES":
                            iline_homo=i_ene_line
                    else:
                        line_ene.append(line)
                        if flag_vline != "YES":
                            iline_homo=i_ene_line
                            #print "alpha homo", iline_homo
                    #print i_ene_line, "here"
                    i_ene_line += 1
        ###### when MO is done, come to here
        ###### 2nd search block :: MO Coeff 
        ###########################################################################
        if MOCoeff_anal==1 and tag_moc=="YES":
            if flag_moc == "NO":
                if re.search(keyword_moc, line) :
                    flag_moc="YES"
            elif flag_moc == "YES" and tag_moc=="YES":
                ## end of MO Coefficients
                if line =='\n':
                    tag_moc="Done"
                    i+=1
                    continue
                #### save MO id
                if imoc==0:
                    #print line, 
                    imo=re.split("[\s\(\)-]+",line)
                    imo=[x for x in imo if x]
                    n_ene_col=len(imo)
                    i+=1
                    imoc+=1
                    continue
                #### save energy
                elif imoc==1:
                    imo_ene=re.split("\s+",line)
                    imo_ene=[x for x in imo_ene if x]
                    del imo_ene[0]
                    #print imo_ene
                    i+=1
                    imoc+=1
                    continue
                ###### arrive MO Coefficients block
                elif imoc <= 2+nbasis:
                    moc_line=re.split("[\s\(\)-]+", line) # coeff's are positive now
                    moc_line=[x for x in moc_line if x]
                    #print moc_line
                    # not storing but display
                    istart = len(moc_line) - n_ene_col 
                    for j in range(n_ene_col):
                        mo_symbol=[]
                        moc_tmp=[]
                        if float(moc_line[istart+j]) > moc_crit:
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
                    ### if read nbasis + 2 line, treat block and initialize
                    if imoc == 1+nbasis:
                        imoc=0
                        imoc_block+=1
                        ### treat block of MOC here
                        #for item in moc_block:
                        #    print item
                        #sorted(moc_block, key=lambda x: x[1])      # this is for tuple
                        #### FORMAT of moc_block: [[P, 1, s], mo_id, mo_ene, mo_coeff
                        sorted_moc_block=sorted(moc_block, key=lambda x:float(x[1]))
                        for item in sorted_moc_block:
                            if print_level >=1:
                                print item
                        k=0     # index for many coeff with the same MO id
                        s='-'
                        id_tmp=0
                        #### SELECT ENERGY LEVEL FOR MPL
                        for eff_coeff in sorted_moc_block:
                            if eff_coeff[0][0] != 'Ni' and eff_coeff[0][0] != 'P':
                                continue
                            if eff_coeff[0][0]=='Ni' and eff_coeff[0][1]=='s':
                                continue
                            #print eff_coeff
                            id_mo=eff_coeff[1]
                            mo_ene=eff_coeff[2]
                            #disp_MOids.append(id_mo)
                            if id_mo == id_tmp:
                                print "\t", s.join(eff_coeff[0]), "\t", eff_coeff[3]
                            else:
                                k=0
                                id_tmp=id_mo
                                #### Save energy for mpl 
                                f_list_ene_select.append(mo_ene)
                                print "\nMO level: ", id_mo, mo_ene
                                print "\t", s.join(eff_coeff[0]), "\t", eff_coeff[3]
                            k+=1
                        #print 
                        moc_block=[]
                        continue
                imoc+=1
                
        ###### 3rd search block :: NBO BonDing extraction
        ##############################################################
        if NBO_anal==1:
            if nfile==3 and fid == 1:
                if flag_nbo == "NO" and tag_nbo == "YES":
                    if re.search(keyword_nbo, line):
                        flag_nbo = "YES"
                ### found NBO and first line                
                elif flag_nbo == "YES" and tag_nbo == "YES":
                    if flag_bd=="NO":
                        if re.search(keyword_nbo_BD, line):
                            flag_bd="YES"
                            bd_line=re.split("[\s\(\)-]+", line)
                            bd_line=[x for x in bd_line if x]
                            print bd_line[0], bd_line[1], bd_line[2], bd_line[3], bd_line[4], bd_line[5], bd_line[6], bd_line[7] 
                            ibd += 1
                        elif re.search(keyword_nbo_end, line):
                            tag_nbo = "Done"
                        else:
                            pass
                    ### in the NBO-BD block, write
                    else:
                        bd_line=re.split("[\s\(\)-]+", line) # - negative sign is eliminated
                        bd_line=[x for x in bd_line if x]
                        if ibd % 2 == 1:
                            print "\t", bd_line[0], bd_line[2], bd_line[3], bd_line[4], bd_line[5], bd_line[6], bd_line[8],
                        else:
                            print "\t",
                            for x in range(len(bd_line)):
                                if float(bd_line[x]) > 0.2:
                                    print  AngSym[x],
                            print
                        if ibd < 4:
                            pass
                        else:
                            flag_bd="NO"
                            ibd=0
                            continue
                        ibd += 1
        i += 1
            
        # end of line: inner line indentation
    # end of readline file
    inf.close() 
    # 
    iline_homo_tmpl.append(iline_homo)
    FL_homo_i_line.append(iline_homo_tmpl)
    if print_level >= 2:
        print list_qchem_ofiles[fid], "homo energy line number: ", FL_homo_i_line[fid]

    FL_ene_line_tmp=[]
    FL_ene_line_tmp.append(line_ene)
    if tag_beta == "YES":
        FL_ene_line_tmp.append(line_ene_b)
    FL_ene_line.append(FL_ene_line_tmp)
    if print_level >= 2:
        print FL_ene_line
    #FL_disp_MO_id.append(disp_MOids)
    FL_ene_select.append(f_list_ene_select)

    fid += 1
###### end of all the files
################################################################################
###### extract energy from MO energy lines
###### using disp_MO ids extract energy
fl_ene_tmp=[]
fid=0
for fline in FL_ene_line: # file index [ f1[line1, line2,...], f2[line1, line2,...
    list_ene=[]
    e_homo_tmpl=[]
    iab=0
    for ab_line_list in fline: # ab index
        i=0        
        #print ab_line_list
    	for e_line in ab_line_list:      # line index, it might be list[0] / [ list[0], list[1]]
    	    ene=re.split("\s+", e_line)
    	    for x in ene:
                if x:       # remove empty entry
	    	        list_ene.append(x)
            if i == FL_homo_i_line[fid][iab]:
                ehomo=list_ene[len(list_ene)-1]
                e_homo_tmpl.append(ehomo)
                #print i, ehomo
            i+=1
        iab += 1
    FL_homo_ene.append(e_homo_tmpl)
    FL_list_ene.append(list_ene)
    fid+=1
print "HOMO energy:: ", FL_homo_ene
if x in  FL_beta_tag:
    print "Orbital energies:: alpha, beta are combined"
'''
# for total MO levels
for f_elist in FL_list_ene:
    print len(f_elist), f_elist
'''
###### Draw mpl
################################ PLOT using MPL ##################################################
### overall variable
fid=0
#print len(FL_list_ene)
### obtain display energy from display list
#for list_ene in FL_list_ene:
for list_ene in FL_ene_select:
    print list_ene
    ### initialize variables for each file drawing
    tag_draw_beta="NO"
    xlength=0.6
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
            ### draw 3 plots for previous energies
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

