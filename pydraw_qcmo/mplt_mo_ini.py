### print option
V_print_link=0
V_print = 0
V_print_mo = 0
V_print_moc = 1
V_print_filter =0
V_print_log = 0
### QChem output file character
Nene_line_QC = 8        # number of energies in a line in Alpha MO block, used in f_extract_ene()

#### key word
key_MOa="Alpha MOs"
key_MOb="Beta MOs"

### key for MO
# first "Occupied" for Alpha, second for beta
KW_MOene_occ="Occupied"
KW_MOene_vir="Virtual"
KW_MOene_beta="Beta MOs"
KW_MOene_end="-----"

### KW for MO Coefficient
KW_MOcoeff="MOLECULAR ORBITAL COEFFICIENTS" # common for alpha & beta
KW_nbasis="basis functions"
Nbasis_show = 2


#### plot option
Nmax_virtual= 5           # normally 3
Nmax_occupy=20
NLINK = 5  

#### MO type for selected atoms
#   ALL     : display all the atoms in the file
#   SEL     : display all the atoms in the list
#   OR      : any atom in SEL
#   AND     : all the atoms in SEL at the same time
#   AND2    : AND 2 groups, first atom and one of the others: 

MO_type_list=['ALL', 'SEL', 'OR', 'AND', 'AND2']

#### NBO
L_NBO_analysis = False

### mplot_qcdraw 
### x divides 3 parts for 3 files
XMIN=0
XMAX=1
# for Ni-CO2 0.05, 1.0 for CO2 -charges
YMAX    =  0.2
YMIN    =  -0.5

"""
if re.search("1-P", files[0]):
    YMAX    =  0.05
    YMIN    = -0.4
elif re.search("2-P", files[0]):
    YMAX    =  0.05
    YMIN    = -0.3
elif re.search("3-P", files[0]):   
    YMAX    =  0.2
    YMIN    = -0.3
else:
    YMAX    = 0.5       # 0.2
    YMIN    = -0.6      # -0.7
"""


