import sys
import re
import os
from atom_valence import *       # fatom_valence
#print "\n".join(sys.path)

######################## MAIN MODULE: draw_MOs_link #####################################
Nmax_lumo=10            # normally 3
Nmax_homo=30
Nlumo_opt=2         # for better display

File_list=[]
#### Can main module give a value to sub_module of draw_MO_ini
#### if not, scan sys.argv again here

#### read qchem output file upto 3
narg=len(sys.argv)
if narg <= 1 :
    print "Error with no arguments: input qchem outfile"
    print "Usage:: sys.argv[0] qchem.out (with job=sp) [L=0|1] [S=a|s]"
    print " in case nfile==3, 1st 3rd are fragments and the middle is complex"
    exit(1)
#### if atom input, selective
#### if not atom input, receive S=option
#### if S=s, obtain atom from filename

for arg in sys.argv:
    if re.search("out", arg):
    	File_list.append(arg)
    if re.search("=", arg):
        option=arg.split('=')
        if re.search("L", option[0].upper()):
            Link_tag=int(option[1])
        elif re.search("S", option[0].upper()):
            Draw_tag=option[1]
        elif re.search("A", option[0].upper()):
            DAtoms=option[1]
        elif re.search("T", option[0].upper()):
            atom_filter_tag=option[1]
if "DAtoms" in locals():
    DAtom_list=re.split("\s",DAtoms)
if 'Draw_tag' not in locals():
    Draw_tag='s'                     # default is 's'elective not 'a'll
if not 'atom_filter_tag' in locals():
    atom_filter_tag='ONE'
#### default for atom for MO drawing    
#else:
#    DAtom_list=["Ni"]
#print DAtom_list
#exit(0)
print "input files in module ", os.path.basename(__file__),  ": ", File_list
Nfiles=len(File_list)

'''
bonding with O
bonding                                             anti-bonding
1                                       88 homo    
2   92(w) 93(w) 94(s,bo,O1) 99(s,bb,C-O1)           100(s,bb-ab,C-O1
3   CO2                                             lumo

homo-1=-1 homo=0, lumo=1, lumo-1=2
MO_link_hl_id=((id_1st_mol, id_2nd_mol),(id_2nd_mol, id_3rd_mol))
'''

#### Fl_nlumo: file_list for maximum num of lumo for each file
#### MO_link_hl_id: pair between 1-2 and 2-3 fragments
if Nfiles == 1:
    Fl_nlumo=[Nmax_lumo]
#### make link if Nfiles >= 2    
elif Nfiles == 2:
    Fl_nlumo=[Nmax_lumo,Nmax_lumo]
    MO_link_hl_id=[]    # does it need for 2 files: CO2 vs bent CO2
elif Nfiles == 3:
    if re.search("1-P", File_list[0]): 
        #MO_link_hl_id=[[[0,1],[0,0],[-1,-1],[-2,-2],[-3,-3],[-4,-4],[-2,-6]], [[1,1],[0,1],[-1,0],[-2,-1],[-3,-1],[-4,0],[-6,-1],[-32,0]]]    # 1-PP w. CO2-bended
        MO_link_hl_id=[[[0,2],[0,0],[-1,-1],[-2,-2],[-3,-3],[-4,-4],[-2,-6]], [[2,1],[0,1],[-1,0],[-2,-1],[-3,-1],[-4,0],[-6,-1]]]    # 1-PP G631s
        #MO_link_hl_id=[[[0,2],[0,0]], [[2,1],[0,1]]]    # 1-PP G631s
        Fl_nlumo=[1,2,2]        # cut number of lumo but include degeneracy of lumo of CO2
    elif re.search("2-P", File_list[0]): 
        MO_link_hl_id=[[[0,0],[0,5]],[[0,1],[5,1]]]
        Fl_nlumo=[2,5,2]        # cut number of lumo but include degeneracy of lumo of CO2
    elif re.search("3-P", File_list[0]):
        MO_link_hl_id=[[[0,0],[0,6]],[[0,1],[6,1]]]
        Fl_nlumo=[2,6,2]        # cut number of lumo but include degeneracy of lumo of CO2
    else:
        Fl_nlumo=[Nmax_lumo,Nmax_lumo,Nmax_lumo]
else:
    print "if num of file is not in between 1 and 3: ", os.path.basename(__file__)
    exit(11)
#print FL_cut_nlumo
#exit(99)

#FL_Atom_types=[             
#    [['Ni','d']],           # 1st type
#    [['C','p'],['O','p']]   # 2nd type
#    ]



#### obtain atom_list then get atom_type in Q-Chee
#### get_atom returns dic or list
#Adic_fname={}
Fl_MO_type=[]
Fl_MO_atoms=[]
Fl_atom_filter_tag=[]

for fname in File_list:
    #### get atom type from filename
    Adic_fname=fatom_valence(fname)
    print "in", os.path.basename(__file__), "atom_dic from filename of", fname,":", Adic_fname

    #### if draw all MO, skip this part
    if re.search('a', Draw_tag.lower()): 
        pass
    #### default: MO_atoms==none, mo_type==none, atom_filter_tag=='ONE', Draw_tag=="'s'
    #### if selective atoms are defined
    elif 'DAtom_list' in locals():
        MO_atoms=DAtom_list
        mo_type=MO_atoms[0]
    #### otherwise obtain atoms from filename or add here in the script
    #### not defined DAtom_list and selectivity

    #### unique for purpose:: three files from input
    #### P for Ni caltalyst -CO2 system and A-fragment
    elif 'P' in Adic_fname.keys() and not 'C' in Adic_fname.keys():
        if re.search("A", fname): ####
    #### P for Ni caltalyst -CO2 system and A-fragment
            if 'N' in Adic_fname.keys():
                #### for 3-PNP
                #mo_type='Ni'
                #MO_atoms['Ni']      # draw for N?
                mo_type='Ni'
                MO_atoms=['Ni']
            else:
                mo_type='Ni'
                MO_atoms=['Ni']      # draw only for Ni
            #MO_atoms=['Ni', 'P']
            #MO_atoms=['N']
            #MO_atoms=['Ni', 'N']      # not working in module
        #### for 1-PP and CO2 complex            
        else:
            mo_type="NiCO"
            MO_atoms=['Ni', 'C', 'O']
            atom_filter_tag="1_OR"     #### Ni is necessary and 1 of C or O
    #### for CO2
    elif 'C' in Adic_fname.keys() and 'O' in Adic_fname.keys():
        mo_type="CO"
        MO_atoms=['C', 'O-1', 'O-2']
        atom_filter_tag="OR"        # if there are 2 of 3 atoms in MO one-id-block
    else:
        mo_type=0
        "Error for Molecular type for selective MO drawing"
    if Draw_tag=="s":
        Fl_MO_type.append(mo_type)
        Fl_MO_atoms.append(MO_atoms)
        Fl_atom_filter_tag.append(atom_filter_tag)
if Draw_tag=="s":
    print "atoms for MO in ", os.path.basename(__file__), Fl_MO_atoms
    print "mo_type in ", os.path.basename(__file__), Fl_MO_type
    print "atom filter tag in ", os.path.basename(__file__) ,": ", Fl_atom_filter_tag
#exit(0)





########################## MODULE: mplot_qcdraw ######################################
#### x divides 3 parts for 3 files
XMIN=0
XMAX=3

# for Ni-CO2 0.05, 1.0 for CO2 -charges
if re.search("1-P", File_list[0]):
    YMAX    =  0.05
    YMIN    = -0.5
elif re.search("2-P", File_list[0]):
    YMAX    =  0.05
    YMIN    = -0.3
elif re.search("3-P", File_list[0]):   
    YMAX    =  0.2
    YMIN    = -0.3
else:
    YMAX    = 0.5
    YMIN    =-2.0

