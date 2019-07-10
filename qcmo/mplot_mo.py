#!/home/joonho/anaconda3/bin/python
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
### 2019. 05. 01 start from dev5, more simplify each part
### 2019. 05. 30 1-PP is done modularized up to BLOCK2
### 2019. 06. 04 degeneracy is covered, beta for Fe is being developed

import re
import os
import argparse
import _math
import numpy as np
from common import *
from mplt_qcdraw import *
from mp_mo_dump import *
from mo_level_link import *
from mplt_mo_ini import *
from qcout_mod import *
from qmo import *
#from mplt_fig_set import *

Xrange_nf = Xrange_nf_fixed_x_length


################### MAIN LOOP ####################
def analyze_files(files,f_atom_lists,sim_atomtype,motypes,base_type,L_link,link_types,link_file,atom_groups,L_NBO,nao_chg_group,nao_chgsum,title,xtitle): 
    print("################################## START Q-CHEM OUTFILE ANALYSIS ##########################")
    print("START {}()::".format(whereami()),files, \
            "\n\t\tMO filter[atom|type] : ", f_atom_lists, motypes, \
            "\n\t\tlink_type :            ", link_types)
    """
        files:          Q-Chem output files, 1, 2, 3, 4, 5 files
        f_atom_lists:   for atom selection in MOC screening
                        if 1: Latom_selection=True
        sim_atomtype:   the same atom but different id such as C|C1 in different files
        motypes:        Choosing MOs using MOC with atomtype such as ONE, SEL, ALL etc
        base_type:      N/A
        L_link:         draw link or not
        link_types:     to calculate link: flow, left|right-split
        link_file:      read link_id.dat and draw link
        atom_groups:    N/A
        title           title in figure
        xtitle          xlabels in figure

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
    Fl_qchem_outf = files[:]        # capital "F" scans all files
    nfile=len(Fl_qchem_outf)
    atomlists_id=expand_dim_str(f_atom_lists)
    ### Link option
    if f_atom_lists:
        Lselect_mo = True
    else:
        Lselect_mo = False
    if link_file:
        if L_link == False: L_link = True

    max_lumo=Nmax_virtual

    ### MO module 
    FL_norbital=[]      # "FL" scans all files with [alpha, beta] components: increasing dimesion of tensor
    FL_homo_ene=[]      # it can have beta, save homo energy [ [f1 a, f1 b],[f2 a, f2 b]... ]
    FL_homo_id=[]       # ene and MO id go parallel
    Fl_homo_ene=[]
    Fl_beta_tag=[]      # 0 for restricted, 1 for beta spin of  unrestricted calculation 

    ### MOC module
    Fl_MO_types=[]      # [ONE(f1), 1sub(f2)...
    Fl_MO_atoms=[]      # just different atoms: [f1: [Ni, P], f2: [C, O], ..
    Fl_MO_atoms_id=[]   # all indices
    Fl_nlumo=[]         # [nlumo-f1, nlumo-f2...    same as a|b
    FL_QCIMO=[]         # list of fl_QCIMO
    FL_MOene_list=[]      # [ [e1, e2... ], [e1, e2...] [...
    FL_MOene_select=[]
    FL_MObcdic_select=[]  # base-coefficient selected
    FL_abMO_dic_list=[]
    FL_MOid_select=[]

    ### LINK part
    FL_sMO_dic=[]
    FL_link_dic=[]
    nhlink =  NLINK

    fid=0
    s='-'       # atom-id linker
    ### Analyze each file
    #####################
    for f in Fl_qchem_outf:
        print (f"############################## FILE {fid+1}/{nfile}:: {f} #########")
        ##### BLOCK 3 :: NBO BonDing extraction #################
        if L_NBO == 1:
            if not isinstance(nao_chg_group[0], list):
                chg_atoms=nao_chg_group[:]
                chgsum_atoms=nao_chgsum[:]
            else:
                chg_atoms=nao_chg_group[fid]
                chgsum_atoms=nao_chgsum[fid]
            with open(f, 'r') as inf:
                analyze_nbo(inf, chg_atoms, chgsum_atoms)
            fid+=1
            if fid == nfile:
                return 0
            else:
                continue
        ### variable block inside file
        #l_QCIMO=[]     # fl_QCIMO=[class list] <- fl_QCIMO=dict{[imo=class QC_IMO]}
        temp_atom_name='Z'      # not atom id such as Ni, Ni-1, Ni-2...
        ##### Variables for Drawing MO Levels
        L_ene_select=[]         
        L_imo_select=[]
        f_ene_dic={}
        f_dic_ab_list=[]
        f_norbital_ab=[]
        disp_MOids=[]

        mo_ene_tmp=-100
        imoc_valuable=[]

        # read one file for analysis & write tmp_file.txt
        lfile = re.split('\.', f)
        tmp0 = 'tmp0-' +lfile[0]+'.log'
        tmp1 = 'tmp1-' +lfile[0]+'.log'
        if fid == 2:
            fname_link = lfile[0]   # save to write link pair

        ####### BLOCK 1 :: MO "Orbital Energies (a.u.)" block
        with open(f, 'r') as inf:
            fL_ene_ab,f_homo_ab,f_homo_id_ab,beta_tag = get_MO_energies(inf)    # len(fL_ene_ab)==2 for the same elevel for a, b
        if vprint_log >= 1: print("----Wrap up MOE")
        ### Wrap up MO Block
        """
            eline_2d =[[Alpha:l1, l2,...],[Beta:l1, l2,...]], FL_homo_ene =[[f1:[A-homo, B-homo],[f2:...
            Fl_homo_ene =[f1-homo, f2-homo,...],              FL_homo_id  =...
        """
        Fl_beta_tag.append(beta_tag)
        FL_homo_ene.append(f_homo_ab)
        if Fl_beta_tag[fid] == 1:
            Fl_homo_ene.append( max(FL_homo_ene[fid][0], FL_homo_ene[fid][1]) )
        else:
            Fl_homo_ene.append( FL_homo_ene[fid][0])
        FL_homo_id.append(f_homo_id_ab)
        FL_MOene_list.append(fL_ene_ab)
        ### keep for energy level only
        f_norbital_ab.append(len(fL_ene_ab[0]))
        if Fl_beta_tag[fid] == 1:
            f_norbital_ab.append(len(fL_ene_ab[1]))
        FL_norbital.append(f_norbital_ab)

        ####### BLOCK 2 :: MOC part (MO Coeff) ########################################
        ### ene & index are obtained here
        ### find key_word_MOcoeff and initialize
        if Lselect_mo :
            with open(f, 'r') as inf:
                L_imo_select,L_ene_select,l_bcdic_select,L_QCIMO,max_lumo,in_mo_atoms,in_mo_atoms_id=obtain_MOC(inf,beta_tag,f_homo_id_ab,f_norbital_ab,atomlists_id[fid],motypes[fid])
            FL_QCIMO.append(L_QCIMO) 
            Fl_nlumo.append(max_lumo)
            if vp_moc >=2: 
                print(f"L_imo_select :: {L_imo_select} after getting from obtain_MOC()")
                print(f"L_ene_select :: {L_ene_select} after getting from obtain_MOC()")


    # for f in files - loop
        ######## FILTER 5: retailing plot list in each file
        #### cut lower energy limit for drawing in figure
        #### make dictionary for each file
        #### change into all levels
        print (f"#################### END: {f}  #######")

        ######### without link option:  just Draw Energy level w/o MOC(link-option)
        #print(f_atom_lists)
        if not Lselect_mo:
            not_digit=[]
            i_lumo_cut= min(f_norbital_ab[0], FL_homo_id[fid][0]+Nmax_virtual)
            i_homo_cut= max(1,FL_homo_id[fid][0]-Nmax_occupy)
            print (f"{whereami()}():: draw all level from {i_homo_cut} {i_lumo_cut}")
            if vprint == 1: print (fL_ene_ab)
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
                L_ene_select.append(tmp_elist[:])
                L_imo_select.append(tmp_list[:])

        if vprint >= 3:
            print ("After imo selection, make dictionary: ", L_imo_select)
            print (L_ene_select)
        print ("<Plot MO levels of ", Fl_qchem_outf[fid], ">")
        print ("   MO_level_ID Energy")
        ### make dictionary
        for imo, moe in zip(L_imo_select[0],L_ene_select[0]) :
            print ("%11d %8.3f" % (imo, moe))
            f_ene_dic[imo]=moe
        f_dic_ab_list.append(f_ene_dic)
        if beta_tag == 1:
            #### change key for alpha, beta mo-id
            #print (L_imo_select[1])
            #print (L_ene_select[1], len(L_ene_select[1]))
            f_ene_dic={}
            for imo, moe in zip(L_imo_select[1],L_ene_select[1]) :
                i_imo=imo
                f_mo=moe
                print (i_imo, moe)
                f_ene_dic[i_imo]=f_mo
            f_dic_ab_list.append(f_ene_dic)
        #### code for beta is developed up to here
        FL_MOid_select.append(L_imo_select)
        FL_MOene_select.append(L_ene_select)
        FL_abMO_dic_list.append(f_dic_ab_list)
        fid += 1
    ######################################################
    ###### end of reading all  the files
    print ("###########################################################")
    print ("############################## Result #####################")
    #print ("       Filename     [HOMO_energy] [ID] beta_tag motype linktype")
    print ("       Filename     [HOMO_energy] [ID] beta_tag linktype")
    print(motypes)
    for i in range(nfile): 
        print ("     %-15s" % Fl_qchem_outf[i], end='' )
        print (FL_homo_ene[i], end='')
        print (" ", FL_homo_id[i], end='')
        print ("%4d" % Fl_beta_tag[i], end='')
        #if motypes:
        #    print ("%8s" % motypes[i], end='')
        if link_types and i > 0:
            print ("%8s" % link_types[i-1])         # NoneType ERROR
        print ("")        # for new line
    print ("################### BLOCK 5: PLOT using MPL ")
    ######################## DRAW MO LEVEL
    #### There might be alpha, beta two dictionaries for one file
    fid=0
    for f_abMO_dic in FL_abMO_dic_list:             # dict = { id: moe, }
        ######## draw one or two dictionary
        beta_tag=Fl_beta_tag[fid]
        if vprint_beta >=1: print(f"Draw beta MO if beta_tag=={beta_tag} for {fid}-th file in {whereami()}()")
        if beta_tag==0:
            # number of files, file id, beta?, alpha or beta, dic[a or b], homo(a or b)
            mplt_ab_draw(nfile, fid, beta_tag, 0, f_abMO_dic[0], float(FL_homo_ene[fid][0]),title=title)
        else:           # 0 for alpha, 1 for beta
            mplt_ab_draw(nfile, fid, beta_tag, 0, f_abMO_dic[0], float(FL_homo_ene[fid][0]),title=title)
            mplt_ab_draw(nfile, fid, beta_tag, 1, f_abMO_dic[1], float(FL_homo_ene[fid][1]),title=title)
            if vp_draw >=2: 
                print(f"{f_abMO_dic[0]} {float(FL_homo_ene[fid][0])} in {whereami()}()")
                print(f"{f_abMO_dic[1]} {float(FL_homo_ene[fid][1])} in {whereami()}()")

        fid+=1
    print("Energy levels for plot is done")
    """
    MO_link_id=(94,99)
                        for left-side link      for right-side link
                     1st-mol 2nd-mol            2nd-mol 3rd-mol
    MOlink_hl_id=[[[ 0       ,0  ],[0 ,   1]], [[0,       1],[    1,     1]]]
    We need          y-val  y-val y-val y-val 
    Link is only working for non-degenerate MOs
    """
    MOlink_hl_id=[[[ 0       ,0  ],[1 ,   1]]]
    MO_link_hl_ene=[]

    ######################### Make MO-link pair
    #################### Obtain link pair from Calculation/link_id file-(made by hand)
    ####### calculate energies for y-axis from imo
    if L_link:
        print ("################### Calculate LINK #####################")
        ######## BASIS-COEFF CHECK
        if vp_link >= 2:
            for i in range(nfile):
                print(f"###### BASE_COEFF's {i+1}:: {Fl_qchem_outf[i]}")
                for qc_class in FL_QCIMO[i][0]:              # jmo, qc
                    attrs = vars(qc_class)              # class:: QC_IMO
                    print (', '.join("%s: %s" % item for item in attrs.items()) )
                    #print(qc_class.basis, qc_class.coeff)
            print("###### BASE_COEFF's END #############")
        ########## LOOP for intervals: nfile-1 interfals for linkage
        link_mo_hl_ids=[]
        n_link_interval=0
        for i in range(nfile-1): 
            ###### f0:f1 in interval
            nlinks = files_nlink(Fl_beta_tag[i], Fl_beta_tag[i+1])
            print(nlinks)
            for f0ab, f1ab in nlinks:
                n_link_interval+=1
                ihomo_f0=FL_homo_id[i][f0ab]       # f0 for left file, f1 for right,
                ihomo_f1=FL_homo_id[i+1][f1ab]     # due to alpha, beta
                #### here obtain link pair
                print (f"############### OBTAIN Link ID in {i+1}-th interval with a-b pairs of {f0ab} {f1ab}" )
                ###### find other MO id's using homo-id
                ###### obtain link pair here
                link_y=[]                       # make link pair
                atomlists_link=atomlists_id[i:i+2]
                #if i < nfile-2:
                #    atomlists_link=atomlists[i:i+2]
                #else:
                #    atomlists_link=atomlists[i:]
            #### OBTAIN link pair from file|Calculation
                if link_file:
                    if n_link_interval==1:
                        fl_molink_hl_id_2D = fread_pair_2Dlist(link_file)
                    fl_molink_hl_id = fl_molink_hl_id_2D[n_link_interval-1]
                    print(f"number of links:: {len(fl_molink_hl_id)}")
                else:
                    if not link_types == None:
                        link_type = link_types[i]
                    else:
                        if nfile == 3 or nfile == 5:
                            link_type, nlink_4atom = get_ltypes_nhlink(nfile, i, atomlists_link)
                        if nlink_4atom:
                            nhlink = nlink_4atom
                    if link_type == 'flow':
                        fl_molink_hl_id=get_link(FL_QCIMO[i][f0ab],ihomo_f0,FL_QCIMO[i+1][f1ab],ihomo_f1,'flow', nhlink, atomlists_link)
                    #### left mo is reactant
                    elif link_type == 'lsplit':
                        #### GET_LINK:: returns link homo-lumo-index
                        fl_molink_hl_id=get_link(FL_QCIMO[i][f0ab],ihomo_f0,FL_QCIMO[i+1][f1ab],ihomo_f1,'lsplit',nhlink,atomlists_link)
                    #### Reverse the L/R side to use right mo, reactant is used as reference
                    elif link_type == 'rsplit':
                        ### input was changed (right mol comes 1st and left mol comes 2nd)
                        atom_rev_lists=[]
                        atom_rev_lists.append(atomlists_link[1])
                        atom_rev_lists.append(atomlists_link[0])
                        fl_molink_hl_id=get_link(FL_QCIMO[i+1][f1ab],ihomo_f1,FL_QCIMO[i][f0ab],ihomo_f0,'rsplit', nhlink, atom_rev_lists)
                    ### write MO Link pair to repair and reuse it
                    link_mo_hl_ids.append(fl_molink_hl_id) 
                    #    f.write fl_molink_hl_id
        # for link_interval 
            # for f0-ab, f1-ab pairs
                ##################### USE ID for draw:: fl_molink_hl_id
                for l_pair in fl_molink_hl_id:
                    if vp_link >= 1:
                        print (l_pair, end='')
                        print("-------------------------")
                    key_f0=ihomo_f0+int(l_pair[0])          # int for reading pairs from file
                    key_f1=ihomo_f1+int(l_pair[1])
                    #print(key_f0, FL_abMO_dic_list[i][f0ab], key_f1, FL_abMO_dic_list[i+1][f1ab])
                    L_pair=True
                    if key_f0 not in FL_abMO_dic_list[i][f0ab]:
                        if vp_draw >=2: print ("Link dictionary error:: exit(44)")
                        #print(FL_abMO_dic_list[i][f0ab])
                        L_pair=False
                        #exit(44)
                    if key_f1 not in FL_abMO_dic_list[i+1][f1ab]:
                        if vp_draw >=2: print ("Link dictionary error:: exit(45)")
                        L_pair=False
                        #exit(45)
                    if L_pair==False:
                        continue
                    #print(f"i={i} f1ab {f1ab} key_f1 {key_f1}")
                    l_ene1=FL_abMO_dic_list[i][f0ab][key_f0]
                    l_ene2=FL_abMO_dic_list[i+1][f1ab][key_f1]
                    link_ene=[l_ene1, l_ene2]
                    link_y.append(link_ene)
                print (" ")
                MO_link_hl_ene.append(link_y)
                #print link_y
    # if L_link: 
        # for i in range(nfile-1): 
            # for f0-ab, f1-ab pairs
        #### when looped all the files using scripts, save link information 
        #if not link_file: fprint_pair_2Dlist("link_id.dat", link_mo_hl_ids)
        #### 3D list, draw with x-coord's
        #if vp_link >= 1:
        #    print ("Link energy:: ")
        #    print (MO_link_hl_ene)
        ##################### DRAW MPL link line
        #print ("Draw Link in ", os.path.basename(__file__))
        #### for each range in MO_link_hl_ene
        #for n in range(nfile-1):
            #nth_links=MO_link_hl_ene[i]
            ### To treat degenerated molecule
                L_degen_f=[]
                for twof in Fl_qchem_outf[i], Fl_qchem_outf[i+1]:
                    fname_list= twof.split('.')
                    if fname_list[0] == 'CO2':  L_degen_f.append(True)
                    else:                       L_degen_f.append(False)
                    
                if not L_degen_f[0]:
                    print("Only no degenracy works for left mol")
                    x=Xrange_nf(nfile, i, Fl_beta_tag[i], f0ab, 1)
                    xl=x[1]                 # 1=right side of 0-th file is left-link
                if L_degen_f[1]:
                    print("Degeneracy 2 works for right mol such as CO2")
                    x=Xrange_nf(nfile,i+1,Fl_beta_tag[i+1],f1ab, 2)
                    xr1, xr2 = x[0], x[2]
                else:
                    x=Xrange_nf(nfile,i+1,Fl_beta_tag[i+1],f1ab, 1)
                    xr=x[0]                             # 0=left side of 1-th file is right-link 
                    
                ### only one set is used here
                x_link=[xl, xr]
                degen=1
                print(f"x-cood's in link: {x_link}")
                #for ene_pair in nth_links:
                for ene_pair in link_y:
                    ### for 5-file double degeneracy
                    #    print("right side has degeneracy")
                    if L_degen_f[0]:
                        pass
                    elif L_degen_f[1]:
                        if degen%2==0:
                            xr=xr1
                            degen+=1
                        else:
                            xr=xr2
                            degen+=1
                    x_link=[xl,xr]
                    if vp_draw >= 1: print(f"{x_link} - {ene_pair} in {whereami()}()")
                    mplt_line_link(x_link, ene_pair)
                            
                    #tmp=ene_pair[1]

    plt.show()
    print ("normal stop"	)


def main():
    parser = argparse.ArgumentParser(description="to analyize QChem outfile and draw MO levels with several files")
    parser.add_argument("-f", "--qcouts", nargs="+", help="more than 1 out files")
    #parser.add_argument("-co", "--mocoeff", action='store_true', help="do coeff analysis")
    #parser.add_argument("-ag", "--atomgroup", nargs='*', help="atom list of 2 groups for mo and link")
    parser.add_argument("-a", "--atomlists", nargs='*', help="atom lists with id for all files")
    parser.add_argument("-at", "--atomtype_sim", nargs='*', help="when the same atoms have different id in two files such as 'C=C1'")
    parser.add_argument("-t", "--motypes", default='ALL', nargs='*', choices=['ONE','ALL','SEL','SUB','OR','AND','AND2'], help="types how to combine the selected atoms also in mo_level_link.py")
    parser.add_argument("-bt", "--basetype", nargs='*', help="add atom basis")
    parser.add_argument("-l", "--link", action='store_true', help="make links for multiple files")
    parser.add_argument("-lt", "--linktypes", nargs='+', choices=['flow', 'lsplit', 'rsplit'], help="for the same molecule (flow) and combination with other molecule (split)")
    parser.add_argument("-lf", "--linkfile", help="read file rather than calculate link")
    parser.add_argument("-sg", "--splitgroup", nargs="+", help="divide atom groups into two group for split link")
    parser.add_argument("-u", "--usage", action='store_true', help="display usage")
    parser.add_argument("--title", help="? will you pass title to mplot lib")
    parser.add_argument("-xt", "--xtitle", nargs='*', help="pass x labels")
    parser.add_argument("--nbo", action="store_true", help="Do NBO analysis")
    parser.add_argument("-ncg","--nao_charge_group", nargs='*', help="for NAO charges and group total charges")
    parser.add_argument("-ncs","--nao_charge_sum", nargs='*', help="for NAO charges sum for a group")
    args=parser.parse_args()

    if args.usage:
        print (" 1 file usage:: %s -f 1-PP-A.out -a \"Ni\" -t SEL" % (os.path.basename(__file__)))
        print (" 2 file usage split:: %s -f 1-PP-A.out 1-PP.out -a \"Ni\" \"Ni C1 O1 O2\" -t SEL SUB -l lsplit" % (os.path.basename(__file__)))
        exit(0)

    #if vprint_log == 1: print ("atomlists:: ", args.atomlists, "in function {}()".format(whereami()))

    analyze_files(args.qcouts,args.atomlists,args.atomtype_sim,args.motypes,args.basetype,args.link,args.linktypes, args.linkfile, args.splitgroup,args.nbo,args.nao_charge_group,args.nao_charge_sum,args.title, args.xtitle)
    return(0)



if __name__ == "__main__":
    main()
