from common import MyClass_str as MyClass
from parsing import str_decom as parse_str

awk                 = MyClass('awk')
backup              = MyClass('backup')
dir_job             = MyClass('dir_job')
git                 = MyClass('git')
start               = MyClass('start')
system              = MyClass('system')
system.install      = MyClass('system.install')
system.manage       = MyClass('system.manage')
system.eta          = MyClass('system.eta')
server              = MyClass('server')
server.mlet         = MyClass('server.mlet')
server.chi          = MyClass('server.chi')
server.kisti        = MyClass('server.kisti')
server.iron         = MyClass('server.iron')
server.ssh          = MyClass('server.ssh')
vmd                 = MyClass('vmd')

start.order     =   "===START Usage===\
                    \n    ORDER:: ~/.bashrc backup ~/bin ...\
                    "
start.bashrc    =   "\n    .bashrc\
                    \n\tturn ON/OFF modules and VM-anaconda\
                    \n\tenvironment modules: intel lammps vasp qchem and so on\
                    \n\tswitch on for anaconda activation\
                    "
start.bin       =   "\n    ~/bin\
                    \n\t$ chi_sys.py\
                    \n\t$ chi_ini.py\
                    \n\t$ research_ini.py\
                    \n\t$ sciwares_ini.py\
                    "
start.menu      =   "\n     comment.py for all system works\
                    \n\t$ show_comment.py\
                    "
start.jobs      =   "\n     comment_subj.py for all jobs\
                    \n\t$ show_comment.py -m subj\
                    "
start.usage     =   "===USAGE of show_comment.py [-h] \
                    \n    $ show_comment.py -m {default:general}\
                    \n    $ show_comment.py -m [general|subj]\
                    \n\tshows all keys\
                    \n    $ show_comment.py -m [general|subj] -j key\
                    \n\tshows details of the key\
                    "
system.start    =   "=== SYSTEM: CHI ==="
system.board    =   "\n    Mother Board::\
                    \n\t# dmidecode | grep 'Product Name'\t#to find out motherboard model"
system.cpu      =   "\n    CPU::\
                    \n\t$ cat /proc/cpuinfo/ \
                    \n\t    Temperature:: sensors \
                    \n\t\t\tinstalled through \"yum install lm_sensors\" \
                    \n\t\t\tconfigured by \"sensors-detect\" \
                    \n\t\t\twhenever check Temperature, command \"sudo sensors-detect\" with \"YES\" for every query, then $sensors\
                    "
system.mem      =   "\n    Memory::\
                    \n\t# dmidecode -t 17 | grep 'Locator\|Size' | grep -v Bank\
                    \n\t    ECC for SERVER\
                    \n\t    4 DDR3 16G equipped\
                    "
system.hdd      =   "\n    HDD::\
                    \n\t# smartctl -a /dev/sd{label}\
                    \n\t# smartctl -t short /dev/sd{label}\twait until calculation done\
                    \n\t# smartctl -l selftest /dev/sd{label}\tcheck the result of selftest\
                    "
system.vga      =   "\n    VGA::\
                    \n\tversion\
                    \n\t    vidia-smi\
                    \n\tmodel\
                    \n\t    lspci –nn | grep VGA\
                    "
system.install.hdd =   "== HDD install ==\
                    \n    Connected: detected in /dev/sd[a|b|c...] \
                    \n    Not Installed: if there is not /dev/sd{label}1, need to make partition    \
                    \n    Install\
                    \n\t     make partition\
                    \n\t\tgdisk /dev/sd{label}\tusing GPT, ext4\
                    \n\t     make filesystem (format)\
                    \n\t\tmkfs -t ext4 /dev/sd{label}\tFormatted\
                    \n\t     manually upload uuid of new partition with mount point in \"/etc/fstab\" \
                    \n\t\tls -al /dev/disk/by-uuid\tThe partitioned hard appear in uuid\
                    \n\t\tvi /etc/fstab\twrite uuid of new partition following the syntax\
                    \n\t\tmount -a\twill read all the mounting points and partitions"
system.install.qchem = "== Q-Chem ==\
                    \n    Installed:: \
                    \n\tv.5 :: /home/joonho/sciwares/archvies/qchem.5.1.tar\
                    "
system.install.intel = "== INTEL ==\
                    \n     Installed:: /opt/intel \
                    \n\tsource /opt/intel/parallel_studio_xe_2019.3.062/bin/psxevars.sh \
                    \n\tCompilation:: \
                    \n\t    serial: $QCHEM_HOME:./configure intel mkl\
                    "
system.manage.user =   "$getent passwd | grep ‘/home’ | cut -d: -f1"
system.eta          = "== ETA ==\
                    \n    HDD:: df -h\
                    \n\tnvm - Mdot2 for debian (250G)\
                    \n\t/sda - Windows10 (250G) not detected\
                    \n\t/sdb /Storage (2.0T): Home_deb, Shared (Data, WCommon)\
                    "

server.chi =         "=== CHI (HOME) ===\
                    \n    Q-Chem::\
                    \n\tsetup .bashrc\
                    \n\t    activate: INTEL\
                    \n\t\tsource ...\
                    \n\t\t\tQ-Chem\
                    \n\t\tsource $QC/bin/...\
                    \n\t$ qchem -np 2 CO2blyp.in CO2blyp.out\
                    \n\t    parallel, but serial run\
                    \n\t$ mpirun -np 4 $QC/exe/qcprog.exe a.in $QCSCRATCH/savename > a.out\
                    \n\t    parallel running with error message\
                    \n\t    check with more number of atoms than nprocess\
                    \n\t    makes 4 scratch folder in $QCSCRATCH\
                    "

server.mlet.system =   "=== MLET (SGE) ===\
                    \n  --System:\
                    \n\tmem 188: node11, 15, 16, 20, 21, 22, 23; 16 nproc * 11.7 G mem\
                    \n\tmem 128: opt07\
                    \n\tmem  96: almost, 12 node * 8G mem\
                    \n\tmem/proc: 2G default, -l mem=nG for qsub\
                    \n\tNode 17, 18, 19 no root passwd, can't run amp\
                    "
server.mlet.plot = "\n  --PLOT Figure\
                    \n\t$ ssh -Y mlet (in login)\
                    \n\t    for drawing in master node\
                    \n\t    not queue-submit job including plot (matplotlib)\
                    \n  --qrun.sh: Running job\
                    \n\te.g.: qrun.sh qname inputfile np nmem\
                    \n\tfix job inside, qchem or vasp\
                    \n\tcall qsub_server.py\
                    \n\t    $ qsub_server.py sge qchem -i input -n np -m mem\
                    "
server.mlet.vasp =   "\n    VASP::\
                    \n\tqsub -N pe500 -pe numa 16 -v np=16 -v dir=pe500 $SB/pypbs/sge_vasp.csh\
                    \n\t    -pe numa: take charge the number of process\
                    \n\tOr Use PBS command\
                    \n\t    qsub_server.py sge \
                    \n\t    qsub_server.py sge -s vasp \
                    \n\t    qsub_server.py sge -s vasp -d dirname -n np[16]\
                    \n\tIN CASE hybrid functional job, it might be killed in 8 hr\
                    \n\t    get node by sleep 'sge.sleep', run at node\
                    "
server.mlet.qchem =  "\n    Q-Chem::\
                    \n\tSetup: .bashrc\
                    \n\t    v5.1\
                    \n\t\tCompilation Method::\
                    \n\t\t    INTEL in server\
                    \n\t\t    MPICH_HOME in /gpfs/opt/openmpi\
                    \n\t\t    $QC=$SCI/qchem5.1p\
                    \n\t\t    $QCAUX = $SCI/qcaux for basis irr. of version \
                    "
server.mlet.sleep   =       "\n    SLEEP::\
                    \n\t$ qsub_server.py sge -s sleep -n 36\
                    \n\t    qsub -pe numa 36 $SB/pypbs/sge_sleep.csh\
                    \n\t$ qsub_server.py sge -s sleep -n 36 -N sleep2\
                    \n\t$ qsub -N amplog -pe numa 36 -l mem=5G $SB/pypbs/sge_sleep.csh\
                    \n\t    : 188G / 36 nproc = 5.22\
                    \n\t$ qsub -N amplog -pe numa 36 -q skylake@node11 $SB/pypbs/sge_sleep.csh\
                    "
server.mlet.at_node =       "\n    RUN @NODE VASP::\
                    \n\t$ sge_vasp_node.csh re0D3mdk_high 36\
                    "
server.sleep        = server.mlet.sleep

server.iron         = "=== Iron:Platinum ===\
                    \n    System::\
                    \n\tlogin node: iron, platinum\
                    \n\thpc cal node: n001 ~ n076\
                    \n\tpartition: X1, X2, X3, X4, X5, n076 (GPU)\
                    \n\tnproc:     8   12  20  24  32\
                    "


server.kisti.pbs    = "=== KISTI ===\
                    \n    System::\
                    \n\tnode(queue)\
                    \n\t    KNL: massive parallel\
                    \n\t\t-l select=4:ncpus=64:mpiprocs=64:ompthreads=1   # max=68, ncpus==mpiprocs\
                    \n\t\t    for VASP NPAR = 8 = sqrt(64)\
                    \n\t\t-q [normal:long]\
                    \n\t\t-l walltime=48:00:00  (for normal)\
                    \n\t\t-l walltime=120:00:00 (for   long)\
                    \n\t    SKL: fast as normal\
                    \n\t\t-l select=1:ncpus=40:mpiprocs=40:ompthreads=1\
                    \n\t\t    for large memory, decrease mpiprocs: ncpus=40:mpiprocs=16\
                    \n\tLocation: PBS scripts\
                    \n\t    Q-Chem:/apps/commercial/test_samples/Q-Chem/qchem_skl.sh, qchem_knl.sh (KNL)\
                    \n    PBS::\
                    \n\tCLI-Options:\
                    \n\t    -N jobname\
                    \n\t    -v fname=fname\
                    \n\t    -v np (N/A)\
                    \n\t    $ qsub -N qjobname -v fname=a[.in] pbs_script.sh\
                    \n\tQ-Chem\
                    \n\t    $ qsub -N CCC6A -v fname=6-CCC-NiFe-A.in qchem_knl.sh\
                    \n\tVASP\
                    \n\t    $ qsub -N dirname $SB/pypbs/pbs_vasp.sh\
                    \n\t\tnumber of process is defined automatically by select*mpiprocs\
                    "
server.ssh.check  = "=== SSH ===\
                    \n    alias.sh::\
                    \n\t$ checknodesamp\
                    \n\t    to check amp_run.sh in running nodes in qstat\
                    \n\t$ checknodes process\
                    \n\t$    scan qstat and check process in the nodes\
                    \n\t$ checknode node process1\
                    \n\t    to check 'python' in running nodes in qstat, defined alias.sh as function\
                    \n\t    $1 can be 'python', 'amp_run.py', 'qchem', etc\n\
                    \n    Scan all the NODES with process name\
                    \n\t$ ssh_mlet_scan_nodes.sh ps process[python]\
                    \n\t$ ssh_mlet_scan_nodes.sh psl process[python]\
                    \n\t    simply denotes the number of the processes\
                    \n    Check One node\
                    \n\t$ checknode 1node process\
                    \n    Kill One node\
                    \n\t$ ssh node pkill process\
                    \n\t$    kills all the processes with the name\
                    \n    Previousely\
                    \n\t$ ssh node01 ps aux | grep process_name(vasp)\
                    \n\t    single node test for vasp\
                    \n\t$ ssh_mlet_scan_nodes.sh process_name[default=vasp]\
                    \n\t    to check (vasp, qcprog, etc) in all nodes\
                    \n    Do process on ONE NODE\
                    \n\t$ ssh_node.sh node_id process_id number_of_processes\
                    \n\t$ ssh_node.sh node13 58412 16\
                    \n\t    echos kill 16 process on node13\
                    \n\t$ ssh_node.sh node13 58412 16 run\
                    \n\t    run kill 16 process on node13\
                    \n    CHECK node for vasp\
                    \n\t$ ssh node08 ps aux | grep vasp | wc -l \
                    "
server.ssh.scripts = "\n    Scripts::\
                    \n\tsee server.pbs.scripts u. -j server -k pbs\
                    "

server.mlet.sge=     "=== SGE: grid-engine (MLET) ===\
                    \n    QSUB\
                    \n\tOverwrite \#$ -option in script\
                    \n\tOptions::\
                    \n\t -N qname (listed in $qstat)\
                    \n\t -pe numa $np (number of process)\
                    \n\t -l mem=10G (memory per process)\
                    \n\t    mem * $np is charged on the node\
                    \n\t -q skylake@node11 (define node), double nodes: -q skylake@node01 -q skylake@node02\
                    \n\t -v vname=value (CLI input can't overwrite the same vname in script)\
                    \n    e.g.\
                    \n\t(Q-Chem) $ qsub -N qname -v qcjob=infile -pe numa np -l mem=3G -v np=np -q skylake@node11 $SB/pypbs/sge_qchem.csh\
                    \n\t(amp)    $ qsub ... double nodes (-q skylake@node01 -q skylake@node02)\
                    \n\t\t\tcores={'node01':ncore,'node02':ncore}\
                    "
server.mlet.qstat=   "\n    SGE command in MLET\
                    \n\tPATH: /gpfs/opt/util\
                    \n\t$ qstat\
                    \n\t    -f: 'qstatf' is aliased\
                    \n\t    -r: details of job\
                    \n\t$ qhist\
                    \n\t$ qfree\
                    \n\t$ qmem\
                    \n\t    available_mem/total_mem(available_proc)\
                    "
server.mlet.sgescripts  = "\n    Scripts::\
                    \n\tSSH\
                    \n\t    ssh_node.sh\
                    \n\t\t...\
                    \n\t    ssh_mlet_scan_nodes.sh\
                    \n\t\tmlet node scan and run (ls rm mkdir vasp qchem ln )\
                    \n\tPBS\
                    \n\t    qstat_ssh.sh\
                    \n\t\tscan qstat s. status r/ then run ps aux/ grep to find execution\
                    \n\t    qrun.sh\
                    \n\t\tto qsub qchem, vasp, amp in MLET\
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
                    \n\t$ ~/bin/backup.sh\
                    \n\t    backup script file\
                    \n\t\t $SB/backup_scripts.py\
                    \n\t    rsync home\
                    \n\t\trsync -avz --delete /home/joonho/ /run/media/joonho/Seagate\ Backup\ Plus\ Drive/Chi_CentOS/home_joonho\
                    \n\t    rsync share_win\
                    \n\t\trsync -avz --delete /shared/share_win/ /run/media/joonho/Seagate\ Backup\ Plus\ Drive/Chi_CentOS/share_win/\
                    \n\t    backup KAIST server: EEWS\
                    \n\t\trsync -avz --delete /Data/EEWS_dat/  /run/media/joonho/Seagate\ Backup\ Plus\ Drive/Chi_CentOS/EEWS_server/\
                    \n\t\trsync -avz --delete /Data/Repository/  /run/media/joonho/Seagate\ Backup\ Plus\ Drive/Chi_CentOS/Repository_chi\
                    "
mountp="/run/user/1000/gvfs/"
backup.camera   =   "\n    BACKUP Phone to Linux Home:\
                    \n\tSamsung galaxy mount point:\
                    \n\t    $mountp\
                    \n\t    mtp is changing whenever connect: /mtp\:host\=%5Busb%3A005%2C004%5D/ \
                    \n\t    $mountp/mtp\:host\=%5Busb%3A005%2C004%5D\
                    \n\t$ rsync -avz /run/user/1000/gvfs/mtp\:host\=%5Busb%3A005%2C004%5D/Card/DCIM/Camera/ /home/joonho/Pictures\
                    "

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
