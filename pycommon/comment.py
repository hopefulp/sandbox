from common import MyClass

water   = MyClass()
vasp    = MyClass()
sge     = MyClass()
pbs     = MyClass()
dir_job = MyClass()
vmd     = MyClass()

water.order = "===WATER===\
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
water.vmdpos2pos = "\n    VMDPOS to POSCAR: vpos_rearrange.py n64.vmdpos -af water_n64.bgf"
water.vasp =        "\n    CONTINUE:: go to 'vasp' attribute"
vasp.order =        "===VASP Usage===\
                    \n    ORDER:: make_incar make_ini run"
vasp.make_incar =   "\n    MAKE incar.key:\
                    \n\t$ vmake_incar -d d3\
                    \n\t\tdefault: -t(dft)=pe, -d(D)=D3\
                    \n\t$ vmake_incar.py -t re0 -d d3\
                    \n\t\thybrid runs with WAVECAR as continous job\
                    \n\t$ vmake_incar.py -t re0 -d d3 -md nvt\
                    \n\t\tto run MD\
                    \n\t$ vmake_incar.py -t revdw\
                    \n\t\tfor revPBE-vdW-DF\
                    \n\t$ vmake_incar.py -t re0vdw\
                    \n\t\tfor revPBE0-vdW-DF    : is this OR?\
                    \n    MAKE INCAR \
                    \n\tvmake_incar.py --rw r\
                    \n\t\tmakes INCAR by reading incar.key with --read option\
                    "
vasp.make_ini =     "\n    MAKE 1st VASP Directory:\
                    \n\t$ vmake_ini.py -a O H -d dirname\
                    \n\t\tKPOINTS=gamma, POTCAR from VaspINI by default and use 'incar.key' for INCAR"
vasp.make_2ndDir =  "\n    MAKE VASP Dir from Dir\
                    \n\t$ vmake_2nd.py pe250s -d pe250 -s POSCAR\
                    \n\t\t make pe250s from dir pe250, default=CONTCAR\
                    "
vasp.run =          "\n    MPIRUN VASP:\
                    \n\t$ mpirun -n 4 ~/sciwares/VASP/vasp.5.4.4/bin/vasp"

sge.vasp =          "===SGE: MLET===\
                    \n    qsub -N pe500 -pe numa 16 -v np=16 -v dir=pe500 $SB/pypbs/sge_vasp.csh\
                    \n\t-pe numa: take charge the number of process\
                    \n    Or Use PBS command\
                    \n\tqsub_server.py sge \
                    \n\tqsub_server.py sge -s vasp \
                    \n\tqsub_server.py sge -s vasp -d dirname -n np[16]\
                    "
pbs.vasp =          "===PBS: KISTI===\
                    \n    qsub -N dirname $SB/pypbs/pbs_vasp.sh\
                    \n\tnumber of process is confirmed in the script 'pbs_vasp.sh'\
                    "
dir_job.clean =     "===DIR: clean, modify filename, jobs for package===\
                    \n    dir_clean_p2.py -j {rm,mv} [-p,-s,-m] -w {qchem, ai, vasp, pbs} -r 'for run'\
                    \n\tjob=rm|mv\
                    \n\twork=qchem|ai|vasp|pbs there will be default for each\
                    \n\t\tpbs to remove *.o\d*, *.e\d*\
                    \n\t-r run for execution\
                    \n\timport common_py\
                    "
dir_job.mod_fname = "\n    dir_fname.py {ls,mvdir,rm,rename} [-p|-s|-m] -a s -i -r\
                    \n\tjob={ls,mvdir,rm,rename}\
                    \n\tmatching [-p|-s|-m]\
                    \n\tappend letter to previous file|dir to use rename\
                    \n\t-i to include directory also\
                    \n\t: rename -p revdw -a s -i\
                    \n\t\tappend -i for all file and directory\
                    "
dir_job.bash =      "\n    CLI_dir.sh\
                    \n\tmodify bash script for simple technique in Command Line Interface\
                    "
vmd.order   =       "=== VMD ===\
                    \n\t:install, general, get single coord(POSCAR), draw box\
                    "
vmd.install =       "\n    INSTALL\
                    \n\tmodify configure\
                    \n\t    :install at HOME\
                    \n\t\t$install_name=\"vmd\"\
                    \n\t\t$install_bin_dir=\"$HOME/.local/bin\"\
                    \n\t\t$install_library_dir=\"$HOME/.local/lib/$install_name\"\
                    \n\t$ ./configure LINUXAMD64\
                    \n\t$ cd src\
                    \n\t$ make install\
                    "
vmd.general =       "\n    General Usage\
                    \n\tFile/New New Molecule...\
                    \n\t    choose where to load file: as New Molecule or on a loaded molecule\
                    \n\t    can select Frames with Stride\
                    \n\t    Browse file and select file type if not automatically determined\
                    "
vmd.save    =       "\n    SAVE one coordinate\
                    \n\t...file/save coordinate\
                    \n\tselect one configuration\
                    \n\tsave as POSCAR format\
                    "
vmd.job_water =     "\n    Water for Lammps and VASP\
                    \n\tLAMMPS: LOAD in 2 step\
                    \n\t    === STEP 1 ===\
                    \n\t\tLOAD BGF\
                    \n\t\t    load a.bgf for configuration\
                    \n\t\t=== Visualization of water\
                    \n\t\t    Graphics/Representations...\
                    \n\t\t\tDrawing Method/DynamicBonds\
                    \n\t\t\t    Bond Radius:0.1, Bond Resolution: 12\
                    \n\t\t=== SET PBC\
                    \n\t\t    vmd> pbc set { a b c alpha beta gamma }\
                    \n\t\t\tget values from bgf file\
                    \n\t\t    vmd> pbc box\
                    \n\t\t\tshow box, if a traj is loaded, this was done at the 1st frame\
                    \n\t\t===Change Display\
                    \n\t\t    Display/Orthographic\
                    \n\t\t===Wrapping: mapping molecule into the box\
                    \n\t\t    vmd> pbc wrap\
                    \n\t    === STEP 2 ===\
                    \n\t\tLOAD Traj(LAMMPS) or OUTCAR(VASP) on the same molecule\
                    \n\t\t    use stride if so many frame\
                    \n\t\t    vmd> pbc wrap -all -compound fragment\
                    \n\t\t\twrap all the frame\
                    "


