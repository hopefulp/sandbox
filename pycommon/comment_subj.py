from common import MyClass
import comment_sys as mod_sys

water   = MyClass()
qcmo    = MyClass()
qchem   = MyClass()
qchem.server = MyClass()
qchem.server.MLET  = MyClass()
qchem.server.CHI   = MyClass()
qchem.server.KISTI = MyClass()
vasp    = MyClass()
vasp.server = MyClass()
vasp.scripts = MyClass()
vasp.postproc = MyClass()
nico2   = MyClass()
myplot  = MyClass()
amp     = MyClass()
amp.server  = MyClass()
amp.scripts = MyClass()
#server.sge     = MyClass()

amp.scripts.run  =  "\n  == Scripts ==\
                    \n\tamp_run.py -f fname -j job -di 1 2 3 4 -nc 4\
                    \n\t    -f : input file=[OUTCAR, extxyz]\
                    \n\t    -j : job=['tr','te','pr'] for training, test, profile; validation was deprecated\
                    \n\t\tprofile draws all PE\
                    \n\t    -di: data index [d1,d2[,d3,d4]] for data region for training and test\
                    \n\t\t-ns: number of data sets are deprecated\
                    \n\t    -ns: training uses n-1/n parts (default n=5)\
                    \n\t         test uses 5th/5 part\
                    \n\t    -nc: Ncore (parallel runs for fingerprints only\
                    \n\t    More Options:\
                    \n\t\t-nc N for number of core for parallel calculation\
                    \n\t\t-tx for use twinx for plot: default=True\
                    \n\t    import myplot2D for plot:\
                    \n\t\tdef draw_dots_two:\
                    \n\t\t   option single y-axis|twinx\
                    \n\t\t   modify color option\
                    \n\t\trefer to myplot\
                    "
amp.scripts.mlet =  "\n\t$ qsub sge_amp.csh\
                    \n\t    qsub script\
                    \n\t    -v scan=ok will run consecutively\
                    \n\t$ sge_amp_scan.py -hl 3 3 3 -nc 12\
                    \n\t    will run scan parallelly\
                    \n\t$ qsub_server.py\
                    \n\t    make qsub command line\
                    \n\t    usage: $qsub_server.py amp -i OUTCAR -qj qname -hl 10 10 -el 0.001 -di 0 1500 1500 2000\
                    \n\t    amp positional argument for software\
                    \n\t    -i input file of PES-force\
                    \n\t    -hl hidden layer, -el energy convergence limit\
                    "
amp.server.chi =   "=== AMP ===\
                    \n  == SERVER ==\
                    \n\tCHI::\
                    \n\t    amp_run.py -f OUTCAR -j tr -nc 4 -di 0 800 1000 -hl 8 8 -el 0.001 -fl 0.01 +g\
                    \n\t\tdetail in amp.scripts\
                    "
amp.server.mlet =   "\n\tMLET::\
                    \n\t    N.B. PLOT ERROR\
                    \n\t\tDo not queue-submit matplotlib\
                    \n\t    PLOT [AMP-Test]:\
                    \n\t\t$ amp_run.py -f OUTCAR -j te -di 1500 2000 -hl 10 10 -el 0.001 -fl 0.1 (-Y master node)\
                    \n\t\t$ amp_plot.py amp_test.txt -f OUTCAR -hl 10 10 -el 0.001 -fl 0.01 -nt 1000 -n 500\
                    \n\t    QSUB [AMP-Training]:\
                    \n\t\t$ qrun.sh amp N1000 OUTCAR 16 4 \"10 10\" 0.001 0.01 \"0 1500 2000\"\
                    \n\t\t    : software_name qname fname np mem=3(G) hl el fl di\
                    \n\t\t    : water N64 mem == 6G\
                    \n\t\t$ qsub_server.py amp -i OUTCAR -qj qname -js te -hl 10 10 -el 0.001 -fl 0.01 -di 0 1500 2000 -m 4\
                    \n\t\t    -qj qname -> qsub -N qname\
                    \n\t\t    -i inputfile\
                    \n\t\t    prints:\
                    \n\t\t    $ qsub -N qname -pe numa 16 -v fname=OUTCAR -v np=16 -v pyjob=tr -v di=\"0 1500 1500 2000\" -v hl=\"10 10\" -v el=0.001 -v fl=0.01 $SB/pypbs/sge_amp.csh\
                    \n\t\t\t -v scan=ok for scan in consecutive way not parallel\
                    \n\t\t:sge_amp.csh\
                    \n\t\t    $PYTHON $EXE -f $fname -j $pyjob -di $di -nc $np -hl $hl -el $el -fl $fl -g\
                    \n\t    submit several jobs with scan\
                    \n\t\t$ sge_amp_scan.py -hl 3 3 3 -nc 12\
                    "
amp.md_anal     =   "\n    MD Analysis\
                    \n\tMD.ene:\
                    \n\t    1st line: \"time    Etot    Epot    Ekin\"\
                    \n\tUsage:\
                    \n\t    $ myplot.py md.ene -x -t MD-Ethylene -yt \"E(eV)\" -xt \"time (10fs)\"\
                    "
                    #\n\t\t    import myplot_default\

myplot.order    =   "===My PLOT===\
                    \n    ORDER:: mplot_start amp_md nico2"
myplot.start    =   "\n    INITIAL SETTING: refer to nico2.mpl_ini\
                    \n\tmyplot_default: moved from myplot2D.py def common_figure():\
                    "
myplot.ini      =   "\n    DEFAULT\
                    \n\tmyplot_default.py\
                    \n\t    figsize: size of figure\
                    \n\t    font.size: rcParams.update()\
                    \n\t    tick_params: labelsize\
                    \n\t    custom_cycler: colors of cyclic order\
                    \n\t    position of text\
                    \n\t\tsingle y-axis: text_x, text_y\
                    \n\t\ttwinx: text_twinx_x, text_twinx_y\
                    "
myplot.md       =   "\n    MD w. AMP\
                    \n\t$ run_plot.py md.ene -x -t MD-Ethylene -yt \"E(eV)\" -xt \"time (10fs)\" \
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
water.vasp_analysis = "\n    ANALYSIS VASP\
                    \n\tgreT a.out\
                    \n\t    step T Etot FreeE Epot Ekin SK SP(?)\
                    \n\tUSE VMD to read OUTCAR\
                    \n\tUSE ASE to read OUTCAR\
                    "
vasp.run            = "=== VASP ===\
                    \n    Structure\
                    \n\t$VASP_HOME/version/bin/vasp, vasp_gam, vasp_std, vasp_ncl\
                    \n\t    : select one following K-points\
                    \n    SERVER\
                    \n\tMLET    \
                    \n\tCHI\
                    \n\t    activate: /opt/intel/compilers_and_libraries_2019.3.199/\
                    \n\t    run: mpirun -np 4 $VASP_HOME/vasp.5.4.4/bin/vasp_std > sidos.out\
                    \n\tKISTI\
                    "
vasp.scripts.make_incar = "\n    === Usage ===\
                    \n\tMAKE incar.key:\
                    \n\t    -sys [bulk|surface|mol]\
                    \n\t    -md [nve|nvt|npt], -t dft[lda,gga,pe,rp,re,re0,revdw,re0vdw,etc]\
                    \n\t    -d dispersion[d2:d3], \
                    \n\t    $ vmake_incar -d d3\
                    \n\t\tdefault: -t(dft)=pe, -d(D)=D3\
                    \n\t    $ vmake_incar.py -t re0 -d d3\
                    \n\t\thybrid runs with WAVECAR as continous job\
                    \n\t    $ vmake_incar.py -t re0 -d d3 -md nve\
                    \n\t\tto run MD\
                    \n\t    $ vmake_incar.py -t revdw\
                    \n\t\tfor revPBE-vdW-DF\
                    \n\t    $ vmake_incar.py -t re0vdw\
                    \n\t\tfor revPBE0-vdW-DF    : is this OR?\
                    \n\tMAKE INCAR \
                    \n\t    vmake_incar.py --rw r\
                    \n\t\tmakes INCAR by reading incar.key with --read option\
                    "
vasp.scripts.make_ini = "\n\tMAKE 1st VASP Directory:\
                    \n\t    $ vmake_ini.py -a O H -d dirname\
                    \n\t\tKPOINTS=gamma, POTCAR from VaspINI by default and use 'incar.key' for INCAR"
vasp.scripts.make_2ndDir =  "\n\tMAKE VASP Dir from Dir\
                    \n\t    $ vmake_d2d.py old_dir new_dir job_type[ini,cont,hybrid,md,dos,band,pchg]\
                    \n\t\t make new_dir from old_dir\
                    \n\t\t ini: copy POSCAR\
                    \n\t\t cont: copy CONTCAR\
                    \n\t\t hybrid: copy WAVECAR etc\
                    " 
vasp.scripts.etc =  "\n\t ase_fconvert.py\
                    \n\t ase_vasp.py\
                    \n\t ase_zpe.py\
                    "
vasp.postproc.p4vasp = "\n    === VASP Post Processig ===\
                    \n\t = P4VASP\
                    \n\t    executable: p4v\
                    "
vasp.postproc.vaspkit = "\n\t = VASP KIT =\
                    \n\t    executable: \
                    "
qchem.server.MLET   = "=== Q-Chem ===\
                    \n    SERVER: \
                    \n\tMLET::   \
                    \n\t    Install:\
                    \n\t\tuse Intel in server\
                    \n\t\t$ ./configure INTEL MKL OPENMPI\
                    \n\t    Setting:\
                    \n\t\t* random /scratch\
                    \n\t\t    qchem script\
                    \n\t\t\tadd QCSCRATCH = `/gpfs/opt/qchem/bin/get_scratch --auto`\
                    \n\t\t* setting mpirun qcparallel.csh\
                    \n\t\tset QCMPI=openmpi in .bashrc\
                    \n\t    Running:\
                    \n\t\t$ qrun.sh qname[dname] fname[.in] nprocess [mem] \
                    \n\t\tfull mem:\
                    \n\t\t    qrun.sh NiFe5    a.in 16 11.7 for full mem 188G\
                    \n\t\t    qrun.sh NiFe5vtz a.in 12 8    for full mem 96G\
                    \n\t\t$ qsub_server.py sge qchem -j CO2M06 -i CO2M06 -n 16 -m 5(G)\
                    \n\t\t$ qsub -N qname -v qcjob=infile -pe numa np -l mem=3G -v np=np -q skylake@node11 $SB/pypbs/sge_qchem.csh\
                    \n\t\tsubmit \"sge_qchem.csh\":\
                    \n\t\t    -v save=ok (in script for qchem save)\
                    \n\t\t    -v iqc=5.1p (in script for qchem version)\
                    \n\t    Check \"howto.sh -j server -k MLET\" for SGE\
                    "
qchem.server.CHI    =   "\n\tCHI::   \
                    \n\t    setup .bashrc\
                    \n\t    $ (parallel) mpirun -np 4 $QC/exe/qcprog.exe a.in $QCSCRATCH/savename > a.out\
                    \n\t\tmakes 4 scratch folder in $QCSCRATCH\
                    "
qchem.server.KISTI  =   "\n\tKISTI-Nurioin::\
                    \n\t$ qsub -N CCC6A -v fname=6-CCC-NiFe-A.in qchem_knl.sh\
                    "
#qchem.server.KISTI_PBS = mod_sys.server.pbs.pbs

qchem.TM        =   "\n    TM: transition metal\
                    \n\tbasis: triple zeta: def2-TZVP, triple zeta valance shell with polarization\
                    \n\tdispersion: Exchange B3LYP\
                    \n\t\t    DFT_D   EMPIRICAL_GRIMME3 \
                    \n\t          : Method wB97X-D3\
                    \n\tUNRESTRICTED    true\
                    "
qchem.out       =   "\n    QCOUT: deal qchem.out file\
                    \n\t$ qcget_georem.pl   -> qcget_in.pl\
                    \n\t$ qcout_geo.pl      a.out : obtain .mol from job=opt\
                    \n\t$ qcget_in.pl       a.out : obtain .mol from job=sp, opt\
                    \n\t$ qcin_mol_rem.pl r=remfile m=molfile [i=outf?] : scratch rem from a.out, mol from b.out(job=sp) make b.in\
                    "


########### PAST WORK ##########################################################################################
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




