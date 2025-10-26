from common import MyClass_str as MyClass
from parsing import str_decom as parse_str
from info_common import filejob
from info_myplot import table
#from common import MyClass
#import comment_sys as mod_sys

amp                 = MyClass('amp')
amp.server              = MyClass('amp.server')
amp.run                 = MyClass('amp.run')
amp.source              = MyClass('amp.source')
amp.scripts             = MyClass('amp.scripts')
amp.alias               = MyClass('amp.alias')
lammps              = MyClass('lammps')
lammps.start            = MyClass('lammps.start')
lammps.lamphet          = MyClass('lammps.lamphet')
lammps.kim              = MyClass('lammps.kim')
myplot              = MyClass('myplot')
packmol             = MyClass('packmol')
qcmo                = MyClass('qcmo')
qchem               = MyClass('qchem')
qchem.server            = MyClass('qchem.server')
qchem.server.mlet       = MyClass('qchem.server.mlet')
qchem.server.chi        = MyClass('qchem.server.chi')
qchem.server.kisti      = MyClass('qchem.server.kisti')
vasp                = MyClass('vasp')
vasp.server             = MyClass('vasp.server')
vasp.scripts            = MyClass('vasp.scripts')
vasp.postproc           = MyClass('vasp.postproc')
vasp.slab               = MyClass('vasp.slab')
vasp.BE                 = MyClass('vasp.BE')
nico2               = MyClass('nico2')
sno2                = MyClass('sno2')
water               = MyClass('water')
mxene               = MyClass('mxene')
hfse2               = MyClass('hfse2')

amp.scripts.run  =  "\n    == Scripts ==\
                    \n\t== AMP direct Run Test\
                    \n\t$ amp_run.py -f OUTCAR -nd 200 -j tr -dt interval -dl 0 150 200 -nc 4 -hl 4 4 -el 0.1 -fl 0.0 +g\
                    \n\t$ amp_run.py -f OUTCAR -nd 1000 -j tr -dt div -dl 5 0 3 -nc 4 -hl 4 4 -el 0.1 -fl 0.0 +g\
                    \n\t$ amp_run.py -f OUTCAR -nd 1000 -j tr -dt div -dl 5 0   -nc 4 -hl 4 4 -el 0.1 -fl 0.0 ! for tr test\
                    \n\t$ amp_run.py -j md -i 0 -ns 50 -f OUTCAR -dt 0.5\
                    \n\t    -f : input file=[OUTCAR, extxyz]\
                    \n\t    -nd: cut ndata from reading input file\
                    \n\t    -j : job=['tr','te','pr','md'] for training, test, profile, MD; validation was deprecated\
                    \n\t\tprofile draws all PE\
                    \n\t    -dt: data selection type ['npart','interval','divide','pick']\
                    \n\t    -dl: data selection list\
                    \n\t\tinterval: [d1,d2[,d3,d4]] for data region for training and test\
                    \n\t\tnpart: number of data sets: training uses n-1/n parts (default n=5) test uses 5th/5 part\
                    \n\t\tdivide: divider, remainer for training, remainer for test\
                    \n\t\tpick: [2, 1] == [ntrain, ntest]\
                    \n\t    -nc: Ncore (parallel runs for fingerprints only\
                    \n\t    MD\
                    \n\t    -i : select the start image from loaded 'OUTCAR'\
                    \n\t    -ns: number of time step\
                    \n\t    -dt: time interval in fs\
                    \n\t    MORE:\
                    \n\t\t-nc N for number of core for parallel calculation\
                    \n\t\t-tx for use twinx for plot: default=True\
                    \n\t    import myplot2D for plot:\
                    \n\t\tdef draw_dots_two:\
                    \n\t\t   option single y-axis|twinx\
                    \n\t\t   modify color option\
                    \n\t\trefer to myplot\
                    "
amp.scripts.mlet =  "\n\t=== AMP Qsub scripts \
                    \n\t$ sge_amp.py -db -i OUTCAR -qj dt15 -m 12G -nc 1 -j tr -hl 8 8 -el 0.001 -fl 0.01 -nt 4766 -ntr 4766 -dt int -dl 0 4766 \
                    \t    : work option\
                    \n\t\t-s for scan\
                    \n\t\t-db for .ampdb\
                    \n\t\tnothing for one-job training in subdirectory\
                    \n\t$ qsub sge_amp.csh\
                    \n\t    qsub script\
                    \n\t    -v scan=ok will run consecutively\
                    \n\t$ qsub_server.py\
                    \n\t    make qsub command line\
                    \n\t    For making all fingerprints\
                    \n\t\t$ qsub_server.py amp -i OUTCAR -qj qname -js tr -hl 10 10 -el 0.001 -dt interval -dl 0 3000 4000 -m 4G\
                    \n\t    amp positional argument for software\
                    \n\t    -i input file of PES-force\
                    \n\t    -hl hidden layer, -el energy convergence limit\
                    "
amp.analysis        = "\n    == Analysis ==\
                    \n\t-- extract Frmse from training results of amp-log.txt\
                    \n\t    grep \"optimization un\" */amp-log.txt -B 1 | awk '{if(NF==11) {print $8} else if(NF==12) {print $9} }'\
                    "
amp.server.chi =   "=== AMP ===\
                    \n    == SERVER ==\
                    \n\tCHI::\
                    \n\t    TR -\
                    \n\t\t$ amp_run.py -f OUTCAR -j tr -hl 8 8 -el 0.001 0.003 -fl 0.01 nt 4000 -ntr 1500 -dtype int -dl 0 3000 +g\
                    \n\t    TEST -\
                    \n\t\trun in ampdb directory\
                    \n\t\t$ amp_run.py -f OUTCAR -j te -tef -hl 5 5 -el 0.001 0.003 -fl 0.0 -nt 4000 -ntr 1500 -dtype int -dl 3000 3002  +g -p EOHL55E0.001F0/amp-untrained-parameters.amp\
                    \n\t\t    -tef: test run force calculation\
                    \n\t\t    detail in amp.scripts\
                    \n\t    MD -\
                    \n\t\t$ amp_run.py -j md -i 0 -ns 50 -f OUTCAR -dt 0.5\
                    "
amp.server.mlet =   "\n\tMLET::\
                    \n\t    N.B. PLOT ERROR\
                    \n\t\tDo not queue-submit matplotlib\
                    \n\t    PLOT [AMP-Test] +g or Unplot at NODE w. -g\
                    \n\t\t$ \"amp_run.py\" shows all\
                    \n\t\t    case 1: amp_run.py -f OUTCAR -j tr -nc 4 -hl 10 10 -el 0.01  -fl 0.1 -nt 5000 -dtype div -dl 5 0 3 +g\
                    \n\t\t    case 2: amp_run.py -f OUTCAR -j tr -nc 4 -hl 10 10 -el 0.001 -fl 0.1 -nt 5000 -dtype int -dl 0 1500 2000 (-Y master node)\
                    \n\t\t$ amp_plot.py amp_test.txt -f OUTCAR -hl 10 10 -el 0.001 -fl 0.01 -nt 1000 -n 500\
                    \n\t    QSUB [AMP-Training]:\
                    \n\t\t1. Make a copy from VASP to amp\
                    \n\t\t    $ make_dir.py dname_new -w amp -j vasp -od vasp_dir\
                    \n\t\t2. Make db\
                    \n\t\t    2.1 Make des-dir\
                    \n\t\t\t$ make_dir.py g2_6 -w amp -j db\
                    \n\t\t\t$ cd g2_6\
                    \n\t\t\t$ sge_amp.py -db -i OUTCAR -qj NN5G2 -nc 10 -j tr -hl 4 -nt 4000 -ntr 100 -dtype int -dl 1000 1100 -m 3G -des gs -pf powNN -pmm 0.05 100.0 -pn 5 -tef\
                    \n\t\t\t$ amp_run.py -f OUTCAR -j tr -hl 4 -el 0.001 -fl 0.01 0.04 -nt 4000 -ntr 100 -dtype int -dl 1000 1100 -des gs -pf powNN -pn 5 -tef \
                    \n\t\t\t    -pf: parameter function: 'log10', 'powNN'\
                    \n\t\t\t    -des: descriptor of 'gaussian', etc\
                    \n\t\t\t    -tef: test for foce which makes files\
                    \n\t\t\tMake amp.db on previous directory by link\
                    \n\t\t\t-m 2G: for database\
                    \n\t\t\t-m 4G: for energy calculation\
                    \n\t\t\t-m 3G: training for 100 image, 12G for 1000 images\
                    \n\t\t\t-nd ndata_total=4000, -dtype data_type=interval -dl data_list=d1~d2 training and d2~d3 test\
                    \n\t\t\tN.B. in nodes, amp_job=tr doesnot runs test; test runs only in master node (login)\
                    \n\t\t    2.11 envs for Module selection\
                    \n\t\t\t(venv) python $SBamp/amp_run.py ......\
                    \n\t\t    2.2 Make db-making directory\
                    \n\t\t\t$ make_dir.py part -w amp -j db\
                    \n\t\t\t$ cd part\
                    \n\t\t\t$ sge_amp.py -db -i OUTCAR -qj Gs26 -nc 10 -j tr -hl 4 -tef -m 3G -LOG10 -DATA \
                    \n\t\t\t$ qsub -N qname -pe numa 16 -l mem=12G -v fname=OUTCAR -v pyjob=tr -v hl='10 10' -v el=0.001 -v fl=0.01 -v ndata=3000 -v dtype=div -v dl='2 0' $SB/pypbs/sge_amp.csh\
                    \n\t\t\t$ amp_run.py -f OUTCAR -j tr -tef -hl 4 -el 0.001 -fl 0.01 -nc 10 -DATA -LOG10\
                    \n\t\t    2.3 Analyze fingerprints\
                    \n\t\t\t$ amp_anal.py -f ../OUTCAR -p amp-untrained-parameters.amp -im 1081 -ia 3 -t 'wrong F'\
                    \n\t\t\t$ amp_anal.py -f ../OUTCAR -p amp-untrained-parameters.amp -im 1081 1083\
                    \n\t\t3.0 Arguments Abbreviation\
                    \n\t\t    -DES:\
                    \n\t\t\tLOG10: -des gs -pf log10 -pmm 0.05 200.0 -pn 10 -pmod del\
                    \n\t\t\tNpowN: -des gs -pf powNN -pn 5\
                    \n\t\t    -DATA:\
                    \n\t\t\t-INT:   -nt 4000 -ntr 1500 -dtype int -dl 1000 [1100]\
                    \n\t\t\t\t-nt 4000 -ntr 100 -dtype int -dl 1000 1100\
                    \n\t\t\t-DIV:   -nt 4000 -ntr 300 -dtype div -dl 3 0\
                    \n\t\t\t\t-nt 3000 -ntr 1500 -dt div -dl 2 0\
                    \n\t\t3. Training\
                    \n\t\t    Run at pwd, Making sub-directory w. single job or scanning\
                    \n\t\t    :single job\
                    \n\t\t\t$ sge_amp.py -qj tr2 -m 3G -nc 10 -hl 4 -el 0.001 -fl 0.1 -tef -DATA -DES  \
                    \n\t\t\t$ amp_run.py -f OUTCAR -j tr -hl 4 4 -el 0.001 -fl 0.1 0.04 -tef -DIV -DES \
                    \n\t\t    :scanning\
                    \n\t\t\t$ sge_amp.py -s -nhl 6 -ihl 1 -qj Fmr -m 12G -nc 8 -hl 5 5 -el 0 -fl 0.01 0.1 -DIV -DES -tef\
                    \n\t\t\t    scanning parameters:\
                    \n\t\t\t    -nhl: number of HL sets\
                    \n\t\t\t    -ihl: interval of number of nodes\
                    \n\t\t4. Test\
                    \n\t\t   test in sub-directory\
                    \n\t\t   run in master node\
                    \n\t\t   $ qrun.sh di te test qname OUTCAR 4 4 '8 8' 0.001 0.00 4000 1500 int '3000 3500' \
                    \n\t\t   $ amp_run.py -f OUTCAR -j te -tef -nc 4 -hl 8 8 -el 0.001 0.003 -fl 0.00 -INT\
                    \n\t\t   $ amp_run.py -f OUTCAR -j te -tef -nc 1 -nt 4000 -ntr 100 -dtype int -dl 1000 1002\
                    \n\t\t   :Scan\
                    \n\t\t   $ sge_amp.py -s -sh 10 -qj EO -j te -tef -m 3G -nc 4 -hl 5 5 -el 0 -fl 0.01 0.1 -DIV\
                    \n\t\t\t    \
                    \n\t\t$ \"qrun.sh\" shows all\
                    \n\t\t    qrun.sh sub_Node $ampjob $sub_dir $qjob     $fin    $np   $mem   \"$hl\"   $el     $fl   ntotal ntrain int \"$data\"\
                    \n\t\t    case 1: qsub training data-interval\
                    \n\t\t\t$ qrun.sh qsub tr  dir   N1500      OUTCAR  16    4  \"10 10\" 0.001   0.01    5000   1500   int  \"0 1500\"\
                    \n\t\t\t  run qsub_server.py\
                    \n\t\t    case 2: direct-run te data-division\
                    \n\t\t\t$ qrun.sh di tr    dir   N1500div32 OUTCAR  16    4  \"10 10\" 0.001   0.00    5000   1500   div  \"3 2\"\
                    \n\t\t\t  run amp_run.py\
                    \n\t\t    : software_name qname fname np mem=3(G) hl el fl di\
                    \n\t\t    : water N64 mem == 6G\
                    \n\t\t$ qsub_server.py amp -i OUTCAR -qj qname -js tr -hl 10 10 -el 0.001 -fl 0.01 -nd 4000 -dt int -dl 0 3000 4000 -m 4\
                    \n\t\t    -qj qname -> qsub -N qname\
                    \n\t\t    -i inputfile\
                    \n\t\t    prints:\
                    \n\t\t    $ qsub -N qname -pe numa 16 -v fname=OUTCAR -v np=16 -v pyjob=tr -v di=\"0 3000 4000\" -v hl=\"10 10\" -v el=0.001 -v fl=0.01 $SB/pypbs/sge_amp.csh\
                    \n\t\t\t -v scan=ok for scan in consecutive way not parallel\
                    \n\t\t:sge_amp.csh\
                    \n\t\t    $PYTHON $EXE -f $fname -j $pyjob -di $di -nc $np -hl $hl -el $el -fl $fl -g\
                    \n\t    submit several jobs with scan\
                    \n\t\t$ sge_amp_scan.py -hl 3 3 3 -nc 12\
                    "                    
amp.run         =   "\n\tRUN::\
                    \n\t    qsub\
                    \n\t\t$ qrun.sh tr pa N2000F1 N2000F1 OUTCAR 16 5 \"10 10\" 0.001 0.1 5000 interval \"0 2000\"\
                    \n\t    ssh node\
                    \n\t\tmake subdir using qrun.sh\
                    \n\t\tcd subdir\
                    \n\t\tamp_run.py -j tr -f OUTCAR -nc 16 -hl 10 10 -el 0.001 -nd 5000 -dt interval -dl 0 1500 2000 -fl 0.1\
                    \n\tTRNINING::\
                    \n\t    qrun.sh\
                    \n\t\tbla\
                    \n\tTEST::\
                    \n\t    qrun.sh\
                    \n\t\t$ qrun.sh di tr N1500divE4 N1500div32 OUTCAR 16 4 '10 10' 0.0001 0.00 4500 1500 div '3 0'\
                    \n\t\t$ qrun.sh di te test N1500 OUTCAR 4 4 '10 10' 0.001 0.00 4500 1500 div '3 0'\
                    \n\t\t$ qrun.sh di te test N1500 OUTCAR 4 4 '10 10' 0.001 0.00 '5000 5455' 1500 int '0 455'\
                    \n\t\t$ qrun.sh di tr N1000divE3F1 N1500div32 OUTCAR 16 4 '10 10' 0.001   0.1  4000  1000   div  '4 0'\
                    \n\t    amp_run.py\
                    \n\t\t$ amp_run.py -f OUTCAR -j te -nt 4500 -ntr 1500 -dt div -dl 3 0 -nc 4 -hl 10 10 -el 0.001 -fl 0.00\
                    \n\t\t$ amp_run.py -f OUTCAR -j te -nt 5000 5455 -ntr 1500 -dt int -dl 0 455 -nc 4 -hl 10 10 -el 0.001 -fl 0.00\
                    \n\t\t$ amp_run.py -f OUTCAR -j tr -nt 4500 -dt div -dl 3 0 -nc 16 -hl 10 10 -el 0.0001 -fl 0.00\
                    \n\t\t$ amp_run.py -f OUTCAR -j tr -nt 4000 -dt div -dl 4 0 -nc 16 -hl 10 10 -el 0.001 -fl 0.1\
                    "
amp.md_anal     =   "\n    == MD Analysis ==\
                    \n\tMD.ene:\
                    \n\t    1st line: \"time    Etot    Epot    Ekin\"\
                    \n\tUsage:\
                    \n\t    $ mplot_f.py -f md.ene -icy2 3 -t AMP-MD -yl \"E(eV)\" -yl2 \"E_kin(eV)\" -xl \"time (10fs)\" -tx \
                    \n\t    $ mplot_f.py -f md_dt01_nt1000.ene -icy2 3 -t AMP-MD -yl 'E(eV)' -yl2 'E_kin(eV)' -xl 'time (0.1fs)' -tx\
                    \n\t\t-icx: x-column, default=1\
                    \n\t\t-tx : twinx\
                    \n\t\t-ity: index for right-y among y-columns\
                    "
amp.source      =   "\n    == ampi ==\
                    \n\t: developed in mlet:~/.local/lib/python3.6/site-packages/amp\
                    \n\t: development in amp.descriptor\
                    \n\t    gs-function parameter generator includes pow(N,sqrt(1/N)), logspace\
                    \n\t    pow(N,sqrt(1/N))\
                    \n\t\tRs is not 0 but moving with eta\
                    \n\t\t    Rc is excluded from parameter generation so Rc is included in source for Rs, eta-rescale: gaussian2.py\
                    \n\t\t\t: control by 'amp_gversion.py'\
                    \n\t\t    Rc is included in parameter generation so Rc is excluded in source: gaussian3.py -- Bug\
                    \n\t\t\t: But cannot be fixed at the moment\
                    "
                    #\n\t\t    import myplot_default\
lammps.start.install ="=== LAMMPS ===\
                    \n    install LAMPHET, KIM, bare-lammps\
                    \n\tLAMPHET compiled by intel_2019: mpirun in $INTEL ...\
                    \n\tKIM     compiled by intel_2019: mpirun in $INTEL ...\
                    "

lammps.lamphet.install ="=== LAMPHET ===\
                    \n   Installation:\
                    \n\tCHI::\
                    \n\t    Intel_2019 [or gcc]\
                    \n\t    after configure, modify Makefile\
                    \n\t\tC[C]FLAGS= ... -fPIC -std=c++11\
                    \n\tMELT::\
                    \n\tKISTI::\
                    "
lammps.lamphet.run  ="\n    RUN::\
                    \n\tPrerequisite\
                    \n\t    lammps input file\
                    \n\t    lammps data file\
                    \n\t    amp.amp - Amp potential\
                    \n\tGenerate PROPhet potential\
                    \n\t    generate_lammps.py\
                    \n\t\tmakes potential_O, potential_H, etc\
                    \n\tRun Lammps\
                    \n\t    Parallel:\
                    \n\t\tmpirun -n 4 $LAMMPS_DIR/lmp_mpi[_shlib] -in in.input -log system.log\
                    \n\tOUTPUT:\
                    \n\t   system.min[heat].lammpstrj, system.log\
                    "
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
                    \n\t$ run_fplot.py -f md.ene -icx 1 -t AMP-MD -yt \"E(eV)\" -xt \"time (0.2 fs)\" \
                    \n\t    -f/-v   : for files or value option\
                    \n\t    -icx    : index of x-value if in file\
                    \n\t    -x      : external input of x-value \
                    \n\t    -tx     : if scale of y-values are different\
                    \n\t    -ity    : for y index for left y-axis\
                    \n\t$ run_fplot.py md.ene -icx 1 -t AMP-MD -yt 'E(eV)' -xt 'time (0.2fs)' -tx -ity 3\
                    "
myplot.nico2    =   "\n    For NiCO2: refer to nico2.myplot"


water.order     =   "===WATER===\
                    \n    ORDER:: calcube makecube pdb2bgf makelmp_in"
water.calcube   =   "\n    CALCULATE CUBE:\
                    \n\t$ chem_math.py -m H2O -d 1.0 -n 64\
                    \n\t    :get cube lattice vector\
                    "
water.makecube =    "\n    MAKE CUBE:\
                    \n\t$ packmol < water_n64.inp\
                    \n\t    :makes a.pdb\
                    \n\tpackmol.inp\
                    \n\t    input(1mol), output with the same type\
                    "
water.xyz2data =    "\n    XYZ to data.lammps\
                    \n\t$ xyz2lmpdata.py water_n768.xyz lattice_a(28.42819)\
                    \n\t    make data.lammps\
                    "
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
packmol.start       ="\n    PACKMOL\
                    \n\tcf. water\
                    \n\t    water.calcube\
                    \n\t    water.makecube\
                    "
packmol.calcube     = water.calcube
packmol.makecube    = water.makecube
packmol.inp         ="\n    XYZ:: Different file format\
                    \n\tTo make xyz, prepare water.xyz file of single molecule\
                    \n\t$ packmol < water_n128.inp\
                    \n\twater_n128.inp::\
                    \n\t    output water_n128.xyz\
                    \n\t    structure water.xyz\
                    "

vasp.run            = "=== VASP ===\
                    \n    Structure\
                    \n\t$VASP_HOME/version/bin/vasp, vasp_gam, vasp_std, vasp_ncl\
                    \n\t    : select one following K-points\
                    \n    SERVER\
                    \n\tMLET    \
                    \n\t    $ qsub -N jobname(queue) -pe numa np -l mem=5G -v np=np -v dir=dirname(log) $SB/pypbs/sge_vasp.csh\
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
vasp.scripts.zpe    = "\n\tContinue to ZPE\
                    \n\t    vas_make_cont.py -p FPt -ex - zpe -j zpe\
                    \n\t\t-p    for prefix make directory\
                    \n\t\t-ex   exlcude for searching '-', 'zpe'\
                    \n\t\t-j    zpe; now works for only zpe\
                    \n\t    Analysis:\
                    \n\t\tzpe_ts_outcar.py dirname -na N ! Only zpe sum for 3*N modes\
                    \n\t    Thermodynamics plot\
                    \n\t\tmplot_table.py FullC54-Ptcp.csv -l -t 'H2 diffusion on Ct-Pt-c'\
                    \n\t\t    import myplot2D --> mplot_levels, mplot_nvector\
                    \n\t\t    input file can deal with a.csv\
                    \n\t\t    -l draws energy level by duplicating the energy\
                    \n\t\t    Usage:\
                    \n\t\t\tmplot_table.py EbPt.csv -xl 'Pt anchoring site'\
                    "
vasp.scripts.etc =  "\n\t ase_fconvert.py\
                    \n\t ase_vasp.py\
                    \n\t ase_zpe.py\
                    "
vasp.BE         =   "\n    == VASP 1st work ==\
                    \n\tmsi: Make structure from Material Studio\
                    \n\tpos2msi.pl {atom series}\
                    \n\t    in case no MD\
                    \n\tpos2msi_mdnew.pl {atom series}\
                    \n\t    in case md colume\
                    "
vasp.postproc.p4vasp = "\n    == VASP Post Processig ==\
                    \n\t = P4VASP + deactivate anaconda to use python2\
                    \n\t    p4v\
                    "
vasp.postproc.vaspkit = "\n\t = VASP KIT\
                    \n\t    vaspkit \
                    "
vasp.server         = "How to run VASP in each server\
                    \n\tUSE: showall.py -j server -k kisti"

sno2.vasp           ="This is for VASP slab band structure\
                    \n\t1. SnO2 slab and surface modification\
                    \n\t2. get POSCAR, cif from material project\
                    \n\t3. make surface using material studio\
                    \n\t4. passivate bottom surface using pseudo hydrogen pp in $VASP_DIR/POTCAR\
                    \n\t5. sort POSCAR with z-axis sort: seperate different H PP\
                    \n\t6. run Opt -> sp (CHGCAR) -> band, dos\
                    \n\t    control sigma of ismear 0 for dos to be more sharp and separatable\
                    \n\t7. Using VBM in bulk\
                    \n\t7.1 Bottom layers passivated w. pseudo-H and 2 bottom layers will be fixed in MD\
                    \n\t    First, plot bottom two layers and find VMB of bottom layers without energy shift\
                    \n\t\tdoslm.py -al 0-23 96-143 24-47 144-191 -ash 72 72\
                    \n\t\t    use l2ldosL2.agr to plot Layer4 & 3 to get absolute energy for band edge\
                    \n\t7.2 Second, plot 4 layer by layer with energy shift\
                    \n\t\tdoslm.py -al 0-23 96-143 24-47 144-191 48-71 192-239 72-95 240-287 -ash 72 72 72 72 -e -1.24\
                    \n\t    cf. doslm.py atom_list atomlist_shape -e VBM -p\
                    \n\t\tatomlist_shape makes group of atomlist\
                    \n\t\t-e shift the reference energy of 0 to VBM\
                    \n\t\t-p for plot\
                    \n\t    Use l2ldos4Lv5_b3.agr for xmgrace format\
                    \n\t    Usage::For sc43\
                    \n\t\tTo find band edge\
                    \n\t\t    $ doslm.py -al 0-23 96-143 24-47 144-191 -ash 72 72\
                    \n\t\tTo draw doscar\
                    \n\t\t    (pristine)\
                    \n\t\t    $ doslm.py -al 0-23 96-143 24-47 144-191 48-71 192-239 72-95 240-287 -ash 72 72 72 72 -e -1.24\
                    \n\t\t    (HOH:to separate OH_B)\
                    \n\t\t    $ doslm.py -al 0-23 96-143 24-47 144-191 48-71 192-239 72-95 240-287 314-317 -ash 72 72 72 70 2 4 -e -1.10\
                    \n\t\t    (ClHNH3:eshift: away -0.92, nearest -0.96, nestN -1.03)\
                    \n\t\t    $ doslm.py -al 0-23 96-143 24-47 144-191 48-71 192-239 72-95 240-287 312-323 -ash 72 72 72 70 2 2 12 -e -0.92\
                    \n\t\t    (HONH4)\
                    \n\t\t    $ doslm.py -al 0-23 96-143 24-47 144-191 48-71 192-239 72-95 240-287 312-325 -ash 72 72 72 70 2 4 10 -e -0.942\
                    \n\t8.3 advanced plot:\
                    \n\t    In case DOS of conduction band is too small\
                    \n\t\tthe DOS in conduction band can be multiplied by fmath.py\
                    "
sno2.fmath          = "pycommon.fmath\n\t" + filejob.fmath 


qchem.server.mlet   = "=== Q-Chem ===\
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
                    \n\t\t$ qsub_server.py qchem -qj Nimono -i CO2M06(.in) -n 16 -m 5(G)\
                    \n\t\t    server is auto-detected\
                    \n\t\t    running one job\
                    \n\t\t$ qsub -N qname -v qcjob=infile -pe numa np -l mem=3G -v np=np -q skylake@node11 $SB/pypbs/sge_qchem.csh\
                    \n\t\tsubmit \"sge_qchem.csh\":\
                    \n\t\t    -v save=ok (in script for qchem save)\
                    \n\t\t    -v iqc=5.1p (in script for qchem version)\
                    \n\t    Check \"howto.sh -j server -k MLET\" for SGE\
                    "
qchem.server.chi    =   "\n\tCHI::   \
                    \n\t    setup .bashrc\
                    \n\t    $ (parallel) mpirun -np 4 $QC/exe/qcprog.exe a.in $QCSCRATCH/savename > a.out\
                    \n\t\tmakes 4 scratch folder in $QCSCRATCH\
                    "
qchem.server.kisti  =   "\n\tKISTI-Nurioin::\
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
                    \n\t$ mplot_f.py -v|-f values|files -j job -t title\
                    \n\t    :: -v y1 y2 y3 ... | -f f1 f2 f3 ...\
                    \n\t    :: -j for job qcmo|ai for xlabel, ylabel, title\
                    \n\t    :: -t, -xt, -yt overwrites xlabel, ylabel, title\
                    \n\t    :: -x for x-column -other options for title\
                    \n\t    --imports my_mplot2d for mplot_nvector\
                    \n\t    --imports plot_job for figure titles for jobs\
                    \n\te.g.:(qcmo) mplot_f.py -v -y -0.8058   -0.7866   -0.7860   -1.0080   -1.2482   -1.2539 -j qcmo -t \"CO2 charges\"\
                    \n\te.g.:(qcmo) mplot_f.py -f nbo-6f.dat -j qcmo -xt Model -ys -1 -yt \"Charge (e)\" -t \"CO2 charges\"\
                    "
nico2.qcmo      =   "\n    QCMO: Plot MO\
                    \n\t$ mplot_f.py -v -y -0.8058   -0.7866   -0.7860   -1.0080   -1.2482   -1.2539 -j qcmo -t \"CO2 charges\"\
                    \n\t$ mplot_f.py -f nbo-6f.dat -j qcmo -xt Model -ys -1 -yt \"Charge (e)\" -t \"CO2 charges\"\
                    "
nico2.eda       =   "\n    EDA: Plot gragh\
                    \n\t$ cd EDA\
                    \n\t$ grep Polar *out | awk '{print $6}'\
                    \n\t$ grep \"CT = DEL\" *out | awk '{print $13}' | tr '\\n' ' '\
                    \n\t$ grep 'SCF Total' *out | awk '{print $13}' | tr '\\n' ' '\
                    \n\t$ cd ..\
                    \n\t$ mplot_f.py -v -y -0.8058   -0.7866   -0.7860  -j eda -t 'CT Energy' -yt 'E (kcal/mol)' \
                    \n\t$ mplot_f.py -f chg-nbo-4f.dat CTene-4f.dat -ys -1 j- -yl 'Charge of CO$_2$ (e$^-$)' 'CT (kcal/mol)' -tx -xv PP PPP PNP PNP-bridged\
                    \n\t$ mplot_f.py -f chg-nbo.dat BE.dat -ys -1 j- -yl 'NAO Charge of CO2 (e$^-$)' 'BE (kcal/mol)' -tx -c r darkcyan\
                    \n\t$ mplot_f.py -f CTene.dat scf.dat -ys 'j-' 'j-' -yl 'CT (kcal/mol)' 'SCF (kcal/mol)' -tx -c red blue\
                    \n\t$ mplot_f.py -f mol5_1froz.dat mol5_2pol.dat mol5_3ct.dat mol5_4scf_total.dat -t EDA -yl 'E (kcal/mol)' -ys j- -yls FRZ POL CT SCF_TOTAL -x Ni5\
                    \n\t$ mplot_f.py -f mol5_nbo.dat mol5_3ct.dat -ys -1 j- -yl 'NAO Charge of CO2 (e$^-$)' -yl2 'CT (kcal/mol)' -tx -x Ni5 -t 'NAO & CT'\
                    \n\t$ mplot_f.py -f mol5_BE-tzD.dat mol5_4scf_total.dat -x Ni5 -ys -1 j- -t 'BE & SCF' -yl 'E (kcal/mol)' -yls BE SCF_TOTAL  \
                    \n\t$ mplot_f.py -f mol5_nbo.dat mol5_BE-tzD.dat -ys -1 -yl 'NAO Charge of CO2 (e$^-$)' -yl2 'BE (kcal/mol)' -tx -t 'BE & NAO' -x Ni5  -yls NAO BE\
                    "
mxene.postjob   =   "Treat series job\
                    \n\tReadme::\
                    \n\t    Functions are kept in alias.sh\
                    \n\t    mv not opt file to a.cont not to be detected\
                    \n\t$ greftail MXNBs22L1*out\
                    \n\t    filename & last energy\
                    \n\t    check filename not to include other files\
                    \n\t$ grefopt MXNBs22L1*out\
                    \n\t    just energy to get energy values\
                    "
mxene.plot      = "Gibbs plot using mplot_gibbs.py"
mxene.myplot    = "\t" + table.mplot_gibbs
mxene.plot2     =  "\tmplot_gibbs.py MXNB-4level.csv -l 'G(U=0)' 'G($U_{Dc}$=1.37)' 'G($U_{Eq}$=2.73)' 'G($U_{Ch}$=3.42)' -c k b g r -ymin -11 -ymax 24\
                   mplot_gibbs.py MX-MXNB-Ueq.csv -l 'MXene' 'MX-NB' -c b r -t 'U$_{Eq}$'\
                   "
hfse2.poscar    = 'POSCAR modification for insertion of new 3-4 O atoms\
                    \n\t$ kpy - is required to run at KISTI\
                    \n\tBombing system temp (500 K) makes the bombing slow -> increase temp by -t\
                    \n\tInsertion at interface needs high T to get over attraction to both sides\
                    \n\t$ kpy pos_modify.py CONTCAR.HfSe2L1O36Hfopt -j add -a O4 -v -t 500 -ht 800 -z 10 -v -d 2 -o HfSe2L1O36HfiO4\
                    \n\t$ kpy pos_modify.py CONTCAR.HfSe2L1O36Hfopt -j add -a O4 -v -t 500 -ht 800 -z 10 -v -vt -z -d 2 -o HfSe2L1O36HfiO4\
                    \nRun VASP\
                    \n\t$ kpy vas_make_ini.py -s POSCAR.HfSe2L1O36MoiO4 -j mdnve -k g -d d2510c\
                    '

def print_obj(job):
    print("Instances:: ", end='')
    if job == None:
        ins_list=[]
        for instance in MyClass.instances:
            ins_name = parse_str(instance)
            if  ins_name not in ins_list:
                ins_list.append(ins_name)
        ### just print the first name
        for ins_name in ins_list:
            print(f"{ins_name}", end=' ')
    print("\n\t    -j for detail")
    return 0


'''
#old version
def print_obj():
    print("Instances:: ", end='')
    for instance in MyClass.instances:
        print(f"{instance}", end=' ')
    print("\n\t    -j for detail")
    return 0
'''




