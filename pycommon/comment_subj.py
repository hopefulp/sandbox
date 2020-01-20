from common import MyClass

water   = MyClass()
qcmo    = MyClass()
nico2   = MyClass()
myplot  = MyClass()
#server.sge     = MyClass()

qcmo.order      =   "===QCMO===\
                    \n    ORDER:: mplot_mo      "
qcmo.mplot_mo   =   "\n    Plot MO using Q-Chem output (job = sp)::\
                    \n\tmplot_mo.py -h\
                    \n\t    call : mplt_ab_draw for drawing\
                    \n\t    Block is modulized in qcout_mod.py\
                    \n\te.g. (level)   : mplot_mo.py -f 1-PP-A.out 1-PP.out 1-PP-B.out \
                    \n\te.g. (link 3file)   : mplot_mo.py -f 1-PP-A.out 1-PP.out 1-PP-B.out -a 'Ni' 'Ni C1 O1 O2' 'C O1 O2' -t ONE SUB ALL -l [-lf m1-3f.dat] \
                    \n\te.g. (link 5file)   : mplot_mo.py -f 1-PP-A.out 1-PP-A.out 1-PP.out 1-PP-B.out CO2.out -a 'Ni' 'Ni' 'Ni C1 O1 O2' 'C O1 O2' 'C O1 O2' -t ONE ONE SUB ALL ALL -l -lf m1-5f.dat  'w. link file'\
                    \n\tParameters::\
                    \n\t\t-l    : draw link, this makes 'link_id.dat', then save to Model1-3f.dat to modify link\
                    \n\t\t-lf   : use link index file, which has 'link1\\nlink2\\nlin...' etc\
                    \n\t\t-a    : atom types to select MO based on the selected atoms_indices -a 'Ni C1 O1 O2' etc\
                    \n\t\t-t --type: how to choose MO based on motype \
                    \n\t\t\tONE if moc-line has one atom in the atom list given by  \
                    \n\t\t\tSEL if moc-line has any of atom \
                    \n\t\t\tALL if moc-line has all the atoms \
                    "
qcmo.mplt_mo_ini=   "\n    Initial condition for plot::\
                    \n\tmplt_mo_ini.py\
                    \n\t    V_print fore verbose\
                    \n\t    YMAX, ... \
                    "
qcmo.qcout_mod  =   "\n    Modularlized Blocks of MO Energies, MO Coefficients, NBO charges \
                    \n\tgeometric average of coefficients of the same base is done at imo_basis_dic()"
qcmo.mo_level_link="\n    Functions::\
                    \n\tget_link(): called by main(); obtain link_id between two files; activated by -l [ltypes] by main\
                    \n\ttweak HOMO LUMO and link"


nico2.order     =   "===Ni-CO2===\
                    \n    ORDER:: mpl qcmo eda"
nico2.mpl_ini   =   "\n    MPL: initialize\
                    \n\t~/.config/matplotlib/matplotlibrc \
                    \n\t:: check by ipython>>>matplotlib.matplotlib_fname()\
                    "
nico2.myplot    =   "\n    Several plots\
                    \n\t$ myplot.py -v|-f values|files -j job -t title\
                    \n\t    :: -v y1 y2 y3 ... | -f f1 f2 f3 ...\
                    \n\t    :: -j for job qcmo|ai for xlabel, ylabel, title\
                    \n\t    :: -t, -xt, -yt overwrites xlabel, ylabel, title\
                    \n\t    :: -x for x-column -other options for title\
                    \n\t    --imports my_mplot2d for mplot_nvector\
                    \n\t    --imports plot_job for figure titles for jobs\
                    \n\te.g.:(qcmo) myplot.py -v -y -0.8058   -0.7866   -0.7860   -1.0080   -1.2482   -1.2539 -j qcmo -t \"CO2 charges\"\
                    \n\te.g.:(qcmo) myplot.py -f nbo-6f.dat -j qcmo -xt Model -ys -1 -yt \"Charge (e)\" -t \"CO2 charges\"\
                    "
nico2.qcmo      =   "\n    QCMO: Plot MO\
                    \n\t$ myplot.py -v -y -0.8058   -0.7866   -0.7860   -1.0080   -1.2482   -1.2539 -j qcmo -t \"CO2 charges\"\
                    \n\t$ myplot.py -f nbo-6f.dat -j qcmo -xt Model -ys -1 -yt \"Charge (e)\" -t \"CO2 charges\"\
                    "
nico2.eda       =   "\n    EDA: Plot gragh\
                    \n\t$ cd EDA\
                    \n\t$ grep Polar *out | awk '{print $6}'\
                    \n\t$ grep \"CT = DEL\" *out | awk '{print $13}' | tr '\\n' ' '\
                    \n\t$ grep 'SCF Total' *out | awk '{print $13}' | tr '\\n' ' '\
                    \n\t$ cd ..\
                    \n\t$ myplot.py -v -y -0.8058   -0.7866   -0.7860  -j eda -t 'CT Energy' -yt 'E (kcal/mol)' \
                    \n\t$ myplot.py -f frozen_1.dat Polar.dat CTene.dat scf.dat -j eda -t EDA -yt 'E (kcal/mol)' -ys -1 -yl FRZ POL CT SCF-TOTAL\
                    \n\t$ myplot.py -f chg-nbo.dat CTene.dat -ys -1 j- -yl 'NAO Charge of CO2 (e$^-$)' 'CT (kcal/mol)' -tx\
                    \n\t$ myplot.py -f chg-nbo-4f.dat CTene-4f.dat -ys -1 j- -yl 'Charge of CO$_2$ (e$^-$)' 'CT (kcal/mol)' -tx -xv PP PPP PNP PNP-bridged\
                    \n\t$ myplot.py -f BE.dat scf.dat -ys -1 j- -t 'BE & SCF' -yt 'E (kcal/mol)' -yl BE SCF-TOTAL\
                    \n\t$ myplot.py -f chg-nbo.dat BE.dat -ys -1 j- -yl 'NAO Charge of CO2 (e$^-$)' 'BE (kcal/mol)' -tx -c r darkcyan\
                    \n\t$ myplot.py -f CTene.dat scf.dat -ys 'j-' 'j-' -yl 'CT (kcal/mol)' 'SCF (kcal/mol)' -tx -c red blue\
                    "
myplot.order    =   "===My PLOT===\
                    \n    ORDER:: mplot_ini amp_md nico2"
myplot.ini      =   "\n    INITIAL SETTING: refer to nico2.mpl_ini\
                    \n\tmy_mplot2d.py:\
                    \n\t    def common_figure():\
                    \n\t\tfigsize: size of figure\
                    \n\t\tfont.size: rcParams.update()\
                    \n\t\ttick_params: labelsize\
                    \n\t\tcustom_cycler: colors of cyclic order\
                    "
myplot.md       =   "\n    MD w. AMP\
                    \n\t$ myplot.py md.ene -x -t MD-Ethylene -yt \"E(eV)\" -xt \"time (10fs)\" \
                    "
myplot.nico2    =   "\n    For NiCO2: refer to nico2.myplot"


water.order =       "===WATER===\
                    \n    ORDER:: calcube makecube pdb2bgf makelmp_in"
water.calcube =     "\n    CALCULATE CUBE:\
                    \n\t$ chem_math.py -m H2O -d 1.0 -n 64"
water.makecube =    "\n    MAKE CUBE:\
                    \n\t$ packmol < water_n64.inp\
                    \n\t\tmakes a.pdb\
                    \n\t\tgopack to see input file"
water.pdb2bgf =     "\n    PDB to BGF:\
                    \n\t$ babel -ipdb water_n64.pdb -obgf water_n64.bgf"
water.modify_bgf = "\n    Modify BGF:\
                    \n\tinclude pbc(CRYSTX); change FF; change charge\
                    \n\tNB: pbc FF coord's bond(CONNECT) are important\
                    \n\tcheck: vmd\
                    "
water.makelmp_in = "\n    Make Lammps Input: \
                    \n\t$ LammpsInput.pl -b water_n64.bgf -f $FF/spcew.ff -s water_n64 -t full\
                    \n\t\tmakes in.water_n64 data.water_n64\
                    \n\t\tcp in.water_n64 data.water_n64 water_n64.bgf to lammps_work_dir"
water.run_lmp =     "\n    RUN Lammps: \
                    \n\t$ mpirun -n 4 ~/.local/bin/lmp -in in.asps -log water.log\
                    "
water.vmd2poscar =  "\n    VMD to POSCAR:\
                    \n\tvmd load a.bgf\
                    \n\tvmd> pbc set { }\
                    \n\tvmd> pbc box\
                    \n\tvmd load b.traj on a.bgf (stride for skip)\
                    \n\tselect one snapshot\
                    \n\tvmd> pbc wrap; move all atoms into the box\
                    \n\tvmd save to poscar\
                    \n\t\t: select only one frame in save panel\
                    "
water.vmdpos2pos =  "\n    VMDPOS to POSCAR: vpos_rearrange.py n64.vmdpos -af water_n64.bgf"
water.vasp      =   "\n    CONTINUE:: go to 'vasp' attribute"
water.vasp_job  =   "\n    VASP JOB for Water\
                    \n\trevPBE+D3 RevPBE0+D3\
                    "
water.vasp_analysis="\n    ANALYSIS VASP\
                    \n\tgreT a.out\
                    \n\t    step T Etot FreeE Epot Ekin SK SP(?)\
                    \n\tUSE VMD to read OUTCAR\
                    \n\tUSE ASE to read OUTCAR\
                    "
