import matplotlib as mpl
import matplotlib.pyplot as plt

################### CONTROL PARAMETERS
### print option
n_loop = 0              # global loop count
vprint = 0
vprint_moe = 0
vp_moc =  1         # up to 3
vprint_log = 0
vprint_filter =0
vprint_loop=0
vp_draw=1
vp_link=0
vprint_beta=0
vp_nbo=0
### QChem output file character
Nene_line_QC = 8        # number of energies in a line in Alpha MO block, used in f_extract_ene()

### KW for MO Coefficient
Nbasis_show = 7        # 2 previously, 5, for Ni-Fe-CO2 7
MOcoeff_crit=0.10                  

#### plot option
Nmax_virtual= 20        # normally 3, 1-PP:5, 2-PPP:8 5-NiFe
Nmax_occupy=100
NLINK = 5               # number of link lines
NH_LINK={'Ni':6, 'O1':2}

#### MPLOT option
mpl.rcParams['xtick.labelsize']=16
mpl.rcParams['ytick.labelsize']=20
mpl.rcParams['lines.linewidth']=5
mpl.rcParams['figure.figsize']=(15, 12)
mpl.rcParams['axes.titlesize']= 'large'
#mpl.rcParams.update({'font.size':30})
XMIN =  0
XMAX =  1

class YvalueModel():
    def __init__(self, ymin3=-0.36,ymin5=-0.40, ymax=0.05):
        self.ymin3=ymin3
        self.ymin5=ymin5
        self.ymax=ymax

Model1 = YvalueModel()
Model2 = YvalueModel()
Model3 = YvalueModel()
Model4 = YvalueModel(ymax=0.20)
Model5 = YvalueModel(ymin3=-0.30)
Model6 = YvalueModel(ymin3=-0.30)

model_ydic={1:Model1, 2:Model2, 3:Model3, 4:Model4, 5:Model5, 6:Model6}
########## Control Y-limit here
Nmodel=5
Nfile=3

YMAX = model_ydic[Nmodel].ymax
if Nfile == 3:
    YMIN = model_ydic[Nmodel].ymin3
elif Nfile == 5:
    YMIN = model_ydic[Nmodel].ymin5

#### number of base for comparision for link: Ni-dxx, O1-px, O2-py etc
NBase_comparison=3


#### MO type for selected atoms
#   ALL     : display all the atoms in the file
#   SEL     : display all the atoms in the list
#   OR      : any atom in SEL
#   AND     : all the atoms in SEL at the same time
#   AND2    : AND 2 groups, first atom and one of the others: 

MO_type_list=['ALL', 'SEL', 'OR', 'AND', 'AND2']

#### NBO
L_NBO_analysis = False


