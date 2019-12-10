from common import MyClass

water   = MyClass()
vasp    = MyClass()
server  = MyClass()
server.sge     = MyClass()
server.pbs     = MyClass()
server.ssh     = MyClass()
backup  = MyClass()
dir_job = MyClass()
vmd     = MyClass()
git     = MyClass()
awk     = MyClass()
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
water.vmdpos2pos =  "\n    VMDPOS to POSCAR: vpos_rearrange.py n64.vmdpos -af water_n64.bgf"
water.vasp =        "\n    CONTINUE:: go to 'vasp' attribute"
water.vasp_analysis="\n    ANALYSIS VASP\
                    \n\tgreT a.out\
                    \n\t    step T Etot FreeE Epot Ekin SK SP(?)\
                    \n\tUSE VMD to read OUTCAR\
                    \n\tUSE ASE to read OUTCAR\
                    "

vasp.order =        "===VASP Usage===\
                    \n    ORDER:: make_incar make_ini run"
vasp.make_incar =   "\n    MAKE incar.key:\
                    \n\t\t-sys [bulk|surface|mol]\
                    \n\t\t-md [nve|nvt\npt], -t dft[lda,gga,pe,rp,re,re0,revdw,re0vdw,etc]\
                    \n\t\t-d dispersion[d2:d3], \
                    \n\t$ vmake_incar -d d3\
                    \n\t\tdefault: -t(dft)=pe, -d(D)=D3\
                    \n\t$ vmake_incar.py -t re0 -d d3\
                    \n\t\thybrid runs with WAVECAR as continous job\
                    \n\t$ vmake_incar.py -t re0 -d d3 -md nve\
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
                    \n\t$ vmake_d2d.py old_dir new_dir job_type[ini,cont,hybrid,md,dos,band,pchg]\
                    \n\t\t make new_dir from old_dir\
                    \n\t\t ini: copy POSCAR\
                    \n\t\t cont: copy CONTCAR\
                    \n\t\t hybrid: copy WAVECAR etc\
                    "
vasp.run =          "\n    MPIRUN VASP:\
                    \n\t$ mpirun -n 4 ~/sciwares/VASP/vasp.5.4.4/bin/vasp"

server.sge.vasp =          "===SGE: MLET===\
                    \n    VASP::\
                    \n\tqsub -N pe500 -pe numa 16 -v np=16 -v dir=pe500 $SB/pypbs/sge_vasp.csh\
                    \n\t    -pe numa: take charge the number of process\
                    \n\tOr Use PBS command\
                    \n\t    qsub_server.py sge \
                    \n\t    qsub_server.py sge -s vasp \
                    \n\t    qsub_server.py sge -s vasp -d dirname -n np[16]\
                    \n\tIN CASE hybrid functional job, it might be killed in 8 hr\
                    \n\t    get node by sleep 'sge.sleep', run at node\
                    "
server.sge.sleep   =       "\n    SLEEP::\
                    \n\t$ qsub_server.py sge -s sleep -n 36\
                    \n\t    qsub -pe numa 36 $SB/pypbs/sge_sleep.csh\
                    \n\t$ qsub_server.py sge -s sleep -n 36 -N sleep2\
                    \n\t    qsub -N sleep2 -pe numa 36 $SB/pypbs/sge_sleep.csh\
                    "
server.sge.at_node =       "\n    RUN @NODE VASP::\
                    \n\t$ sge_vasp_node.csh re0D3mdk_high 36\
                    "
server.ssh.nodes   =       "=== SSH ===\
                    \n    Scan all the NODES for process name\
                    \n\t$ ssh node01 ps aux | grep process_name(vasp)\
                    \n\t    single node test for vasp\
                    \n\t$ ssh_sge_nodes.sh process_name[default=vasp]\
                    \n\t    to check (vasp, qcprog, etc) in all nodes\
                    "
server.ssh.node    =       "\n    Do process on ONE NODE\
                    \n\t$ ssh_node.sh node_id process_id number_of_processes\
                    \n\t$ ssh_node.sh node13 58412 16\
                    \n\t    echos kill 16 process on node13\
                    \n\t$ ssh_node.sh node13 58412 16 run\
                    \n\t    run kill 16 process on node13\
                    "
server.ssh.check_nod =     "\n    CHECK node for vasp\
                    \n\t$ ssh node08 ps aux | grep vasp | wc -l \
                    "

server.pbs.vasp =          "===PBS: KISTI===\
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

git.order       =   "=== GIT ===\
                    \n\t:push pull pull w. force remote \
                    "
git.push        =   "\n    PUSH\
                    \n\t$ git add . -A\
                    \n\t$ git commit -m \"message\"\
                    \n\t$ git push [origin master]\
                    "
git.pull        =   "\n    PULL\
                    \n\t$ git pull\
                    \n\tIf not working, use forced pull\
                    "
git.overwrite   =   "\n    FETCH & RESET\
                    \n\toverwrite the local changes: if believe the stage in HEAD\
                    \n\t$ git fetch --all\
                    \n\t$ git reset --hard origin/master\
                    "
awk.vasp_logfile =  "=== AWK ===\
                    \n    more job.remdk | awk '{ if($4) {split($4,arr,\":\"); print arr[1]-hour, arr[2]-min;} {hour=arr[1]; min=arr[2];}}'\
                    "
backup.crontab  =   "=== BACKUP ===\
                    \n    CRONTAB to backup hd: \
                    \n\toptions: \
                    \n\t    :: -e for edit\
                    \n\t    :: -l for list\
                    \n\te.g.\
                    \n\t    ::min hr dat mon day command\
                    \n\t    :: 0  1   *   *   *   rsync -avz --delete /home/joonho/ /NAS1/home_joonho\
                    \n\t    :: 0  5   1   *   *   rsync -avz --delete /NAS1/home_joonho/ /NAS2/bak_home/\
                    "
onedrive="/c/Users/hopef/OneDrive"
remote_dir="joonho@chi.kaist.ac.kr:/Shared_win/shared_win/Windows10_backup"
backup.windows  =   f"\n    Windows10 BACKUP to Linux\
                    \n\tInstall::\
                    \n\t    Git Bash\
                    \n\t    rsync\
                    \n\t\tuse -avzz for compression option\
                    \n\tWORK::\
                    \n\t    open Git Bash\
                    \n\t    $ cd /c/WinData\
                    \n\t    $ ./rsync_win.sh\
                    \n\trsync error\
                    \n\t    protocol version mismatch -- is your shell clean?\
                    \n\t\t\"turn off anaconda in server chi\" \
                    \n\trsync_win.sh\
                    \n\t    rsync -avzz --delete /c/WinData/            {remote_dir}/WinData\
                    \n\t    rsync -avzz --delete {onedrive}/Documents/   {remote_dir}/win_docu\
                    \n\t    rsync -avzz --delete {onedrive}/Pictures/    {remote_dir}/win_pic\
                    "
backup.hard     =   "\n    BACKUP to External hd:\
                    \n\t    backup script file\
                    \n\t\t/home/joonho/sandbox_gl/backup_scripts.py\
                    \n\t    rsync home\
                    \n\t\trsync -avz --delete /home/joonho/ /run/media/joonho/Seagate\ Backup\ Plus\ Drive/Chi_CentOS/home_joonho\
                    \n\t    rsync share_win\
                    \n\t\trsync -avz --delete /shared/share_win/ /run/media/joonho/Seagate\ Backup\ Plus\ Drive/Chi_CentOS/share_win/\
                    \n\t    backup KAIST server: EEWS\
                    \n\t\trsync -avz --delete /Data/EEWS_dat/  /run/media/joonho/Seagate\ Backup\ Plus\ Drive/Chi_CentOS/EEWS_server/\
                    \n\t\trsync -avz --delete /Data/Repository/  /run/media/joonho/Seagate\ Backup\ Plus\ Drive/Chi_CentOS/Repository_chi\
                    "
                    
