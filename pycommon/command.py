#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_files, MyClass, list2str
import comment_subj
import comment_sys
from server_env import nXn
from datetime import datetime

amp     =   MyClass('amp')
plot    =   MyClass('plot')
ga      =   MyClass('ga')
nc      =   MyClass('nc')
dac     =   MyClass('dac')
slurm   =   MyClass('slurm')
kisti   =   MyClass('kisti')
pbs     =   MyClass('pbs')

user=os.getenv('USER')


### modify input-detail here
### AMP-NODE job
amp_data_type={'interval': '-nt 4000 -ntr 100 -dtype int -dl 1000 1100 1200',
                'interval-te': '-nt 4000 -ntr 100 -dtype int -dl 1000 1100',
                'division': '-nt 4000 -ntr 1000 -dtype div -dl 4 0'}
amp_gaussian_function_param={'log': ' -des gs -pf log10 -pmm 0.05 200 -pn 10 -pmod del', 'pow': ' -des gs -pf powNN -pn 5'}

amp_node_str={ 'noforce': ' -j tr -hl 4 -el 0.001 -fl -0.1  -nc 10',
                'force' : ' -j tr -hl 4 -el 0.001 -fl 0.1 0.04 -nc 10'
                }
amp_node_str_suff=' -hf -g'     # hf: hack force, -g: suppress drawing graph
### AMP-SERVER job
amp_qsub_str={ 'force' : "-pe numa 10 -l mem=3G -v fname=OUTCAR -v pyjob=tr -v hl=4 -v el=0.001 -v tef=force -v fl='0.1 0.04'",
               'noforce': "-pe numa 10 -l mem=3G -v fname=OUTCAR -v pyjob=tr -v hl=4 -v el=0.001 -v tef=force -v fl='-0.1' "
               }
amp_qsub_data_type={'interval':  "-v nt=4000 -v ntr=100 -v dtype=int -v dlist='1000 1100 1200'",
                    'division': "-v nt=4000 -v ntr=1000 -v dtype=div -v dlist='3 0'"
                    }
amp_qsub_gaussian_param={'log': "-v des=gs -v pf=log10 -v pmod=del -v pmm='0.05 200.0' -v pn=10",
                        'pow': "-v des=gs -v pf=powNN -v pn=5"
                        }

### AMP-GET-QSUB job
amp_getqsub_str = { 'noforce': ' -j tr -hl 4 -el 0.001 -fl -0.1', 'force': ' -j tr -hl 4 -el 0.001 -fl 0.1 0.04 -nc 10'}
amp_getqsub_datatype = amp_data_type
amp_getqsub_symmfuncttype = amp_gaussian_function_param


amp.clean       =   "\nCLEAN\
                    \n    clean.py -w amp pbs -a -y\
                    "
amp.statistics  =   "\nSTATISTICS\
                    \n    stat_check.py -st [e|f|ef]\
                    "
amp.ga          =   "\nGA: AMP\
                    \n    ga_amp_run.py -js qsub -nch 12 -hl 7 -nn 15 -np 4\
                    "

amp.analysis    =   "\nAnalysis\
                    \n    Training set: amp-log.txt\
                    \n\tget_TRFrmse NN 2 \
                    \n\tget_TRFrmse NN subdir 2\
                    \n\tget_TRFrmse args in $pyamp/amp_bashfunc.sh\
                    \n\t    null\t: scan all the directory with amp-log.txt and get training results\
                    \n\t    1\t:with one column of value\
                    \n\t    NN[dirname]: prefix or dirname\
                    \n\t    NN 2\t: including directory name\
                    \n\t    use with sort for plot\
                    \n    Test set: test_fstat_acc.txt\
                    \n\tget_TEFrmse NN subdir 2\
                    "
ga.amp          =   amp.ga
ga.dyn          =   "\nGA: Dynamics\
                    \n    ga_dyn_run.py -e 'mllorenz.py' -js node -mi 200000 -nd 160000 -hl 8 -nn 15\
                    \n\tpytorch+gpu coding\
                    \n    GPU: iron, n076\
                    "

pbs.queue       =   f"\tqs = qstat -u {user} ! alias\
                    \n\th5 = watch -d -n 300 'ls -alt' ! to protect the connection in cf. KISTI\
                    \n\thk = watch -d -n 60 'kpy qst.py'\
                    "

slurm.pbs       =   pbs.queue +  "\n\tpe     to get free processes\
                    \n\tpef    to get free nodes\
                    \n\tpestat to see all nodes\
                    \n\tsq: to check pending jobs in each partition\
                    "
slurm.amp       =   " AMP in slurm\
                    \n\t1. amp_jobs.sh fp outcar npartition idata | sh\
                    \n\t    amp_mkdir.py $dname -w amp -j des\
                    \n\t\tprepare pure directory by coping OUTCAR\
                    \n\t\tcalculate 1st data section resulting ampdb's\
                    \n\t1.9 amp_jobs.sh db OUTCARK1000Kdt.5 2 60 | sh\
                    \n\t    Make feature directory such as NN10: serial job in master node\
                    \n\t\tlook up to the 'amp_jobs.sh' to find out the list of features\
                    \n\t\tamp_mkdir.py dirname -w amp -j db\
                    \n\t\t    which makes NNfeature directory\
                    \n\t\tamp_wrapper.py -f $outcar -js sbatch -j tr -qn $dname -p $npart -sf NN$num -dl $ndata  $(expr $ndata + 120) -t db &\
                    \n\t\tsubmit: sbatch -J NN5tr -p X2 -N 1 -n 12 sbatch_tr.csh\
                    \n\t\t    check wrapper code to submit to slurm\
                    \n\t2. amp_jobs.sh db OUTCARK1000Kdt.5 1 | sh\
                    \n\t    Make full db as for all the OUTCAR's\
                    \n\t\tdb outcar npartition\
                    \n\t\tdata list should be given in the script\
                    \n\t3 Training\
                    \n\t    amp_jobs.sh wrapper2 OUTCAR 2 | sh\
                    \n\t\tset dirname symmetry-function inside script\
                    \n\t    comment out idipol in VASP cal for MD, which makes energy pump -> extract positive energy in data preprocessing\
                    "


plot.first      =   "\nPLOT\
                    \n    TRAIN & TEST\
                    \n\tDraw file\
                    \n\t    f2d_plot.py train.dat -pd -xs 2\
                    \n\t    f2d_plot.py train.dat test.dat -pd -t 'TR-100:TE-100'\
                    \n\t\t-pd pandas, -xs xspacing, -xst xspacing type (default:numerically)\
                    \n\tamp_jobs.sh grep\
                    \n\tampplot_dir.py -p NN -t 'Ndata 100 (Gs-pow)'\
                    \n\tampplot_dir.py -p NN -t 'Training Set: Gs-pow' -y tr -yd . Nd300hl2020 Nd500hl2020\
                    \n\tampplot_dir.py -p NN -t 'Test Set    : Gs-pow' -y te -yd . Nd300hl2020 Nd500hl2020\
                    "
nc.build        =   "    BUILD\
                    \n\tnc_build.py build -n graphene -sn rect -sc n m -g -s -o name -fm\
                    \n\t    build, etc\
                    \n\t    -n name [graphene|]\
                    \n\t    -sn subname [hexa,rect,gnrac,gnrzz]\
                    \n\t    -sc size of supercell\
                    \n\t    -g view structure via xcrysden\
                    \n\t    -s save structure\
                    \n\t    -o outfile\
                    \n\t    -fm format [vasp(POSCAR)|\
                    "
dac.build       =   "    Build graphene using nanocore\
                    \n\tSort poscar\
                    \n\t    pos_sort.py POSCAR.NiNi -mv\
                    \n\t\t(loc) pyvasp\
                    \n\t\t-mv write to the input file\
                    "

kisti.shell     =   "\t:: Commands in alias.sh or PYTHONPATH\
                    \n\t(kpy) to run python\
                    \n\t    kpy to run default python instead of shebang\
                    \n\t    $kpy command.py kisti\
                    \n\t\t: the same as 'python $(which command.py ) kisti\
                    "
kisti.pbs       =   f"\tqst.py : runs 'qstat -f' to see long jobnames\
                    \n\t    $ kpy qst.py\
                    \n{pbs.queue}\
                    "

def show_command(work, subwork, job_submit, qname, vjob, vsjob, inf, keyvalues, nodename, nnode, nproc, nodelist, sftype, dtype, partition,poscar, nhl,idata,ndata):
    
    if keyvalues:
        kw1 = keyvalues[0]
    ###### make command line in advance
    ### KISTI
    kisti.vas = f"\t(VASP) :: (std, skylake, pbs_vasp.sh --> pbs_vasp_kisti_skl.sh)\
                \n\t\t$ kpy vas_make_ini.py -s POSCAR.{qname} -al : all prepared except POSCAR\
                \n\t\t$ kpy vas_make_ini.py -s POSCAR.{qname} -j {vjob}\
                \n\t\t$ qsub -N {qname} $SB/pypbs/pbs_vasp.sh "
    kisti.vas += f"\n\t    :: RERUN Opt for failed opt"
    kisti.vas += f"\n\t\t$ qsub -N {qname} $SB/pypbs/pbs_vasp_kisti_sklopt.sh"
    kisti.vas += f"\n\t    :: KPOINTS Sampling"
    kisti.vas += f"\n\t\t$ python $(which vas_make_ini.py) -j kp -s POSCAR -kd 2 -kps 2 1 4 1 6 1 8 1"
    kisti.vas += f"\n\t    :: FAST Run using Backfillingn"
    kisti.vas += f"\n\t\t$ qsub -N {qname} -l walltime=1:00:00 $SB/pypbs/pbs_vasp.sh "  
    kisti.vas += f"\n\t    :: SLAB XY-relax"
    kisti.vas += f"\n\t\t$ kpy vas_make_ini.py -s POSCAR.{qname} -j copt -exe xyrelax"
    kisti.vas += f"\n\t\t$ qsub -N {qname} -v crelax=yes $SB/pypbs/pbs_vasp_kisti_skl.sh"
    kisti.vas += f"\n\t    :: long"
    kisti.vas += f"\n\t\t$ qsub -N {qname} $SB/pypbs/pbs_vasp_kisti_skllong.sh"
    kisti.vas += f"\n\t\t$ qsub -N {qname} -q long -l walltime=96:00:00 $SB/pypbs/pbs_vasp_kisti_skl.sh"
    kisti.vas += f"\n\t    :: GAMMA, ncl"
    kisti.vas += f"\n\t\t$ qsub -N {qname} -v exe=gamma $SB/pypbs/pbs_vasp_kisti_skl.sh"
    kisti.vas += f"\n\t    :: MEMORY Issue run half process to save memory usage"
    kisti.vas += f"\n\t\t$ qsub -N {qname} -l select=20:ncpus=40:mpiprocs=20:ompthreads=1 $SB/pypbs/pbs_vasp_kisti_skl.sh"
    kisti.vas += f"\n\t\t$ qsub -N {qname} -l select={nnode}:ncpus=40:mpiprocs={nproc}:ompthreads=1 $SB/pypbs/pbs_vasp_kisti_skl.sh"
    kisti.vas += f"\n\t\t$ qsub -N {qname} $SB/pypbs/pbs_vasp_kisti_skl2.sh"
    kisti.vas += f"\n\t\t    pbs_vasp_kisti_skl2 for half use of cpu for memory issue"
    kisti.vas += f"\n\t    :: FAKER Job & OVERwrite"
    kisti.vas += f"\n\t\t$ kpy vas_make_ini.py -j fake -s POSCAR.{qname} -sj {vjob} -al -ra -d d{datetime.now().strftime('%d%H')} -n 6 -e g : more info_vasp.py"
    kisti.vas += f"\n\t\t$ kpy vas_make_ini.py -j {vjob} -r on -s POSCAR.{qname} -d dnameid : o-overwrite n-not submit job"
    kisti.vas += f"\n\t    :: (OORINano)"
    kisti.vas += f"\n\t\t$ qsub -N {qname} -v cat='orr' pbs_vasp_kisti_skl.sh"
    ### IRON(slurm)
    if not nproc:
        nproc = nnode * nXn[partition]
    if poscar and re.match('POSCAR', poscar) :
        if poscar[7:] != 'name':
            dirname = poscar[7:]
    if 'dirname' not in locals():
        dirname = qname
    ncpu =  int(nXn[partition]/2)
    slurm.vas = f"==== Job submission\
            \n\t:: INI\
            \n\t    $ vas_make_ini.py -s POSCAR.{dirname} -j {vjob} -x {partition} -N {nnode} -np {nproc}\
            \n\t    $ vas_make_ini.py -s POSCAR.{dirname} -j {vjob} -x {partition} -N {nnode} -np {nproc} -al : all prepared except POSCAR\
            \n\t    $ vas_make_ini.py -s POSCAR.{dirname} -j {vjob} -sj {vsjob} -x {partition} -N {nnode} -np {nproc} -al : all prepared except POSCAR\
            \n\t    : if qsub fails but dir exists\
            \n\t\tsbatch -J {dirname} -p X{partition} -N {nnode} -n {nproc} /home/joonho/sandbox/pypbs/slurm_sbatch.sh\
            \n\t\tsbatch -J {dirname} -p X{partition} -N {nnode} -n {nproc} --export=exe='gam' /home/joonho/sandbox/pypbs/slurm_sbatch.sh\
            \n\t\tsbatch -J {dirname} -p X{partition} -N {nnode} -n {nproc} /home/joonho/sandbox/pypbs/slurm_sbatch_sim.sh\
            \n\t    : to run vasp.5\
            \n\t\tsbatch -J {dirname} -p X{partition} -N {nnode} -n {nproc} --export=v=5 /home/joonho/sandbox/pypbs/slurm_sbatch.sh\
            \n\t:: kpoints sampling\
            \n\t    $ vas_make_ini.py -j kp -s POSCAR -kd 2 -kps 2 1 4 1 6 1 8 1 -x {partition} -N {nnode}\
            \n\t:: For memory issue\
            \n\t\tsbatch -J {dirname} -p X{partition} -N {nnode} -c {ncpu} --export=hmem=1 /home/joonho/sandbox/pypbs/slurm_sbatch.sh\
            \n\t\tsbatch -J {dirname} -p X{partition} -N {nnode} --ntasks-per-node {ncpu} --export=hmem=1 /home/joonho/sandbox/pypbs/slurm_sbatch.sh\
            \n\t:: continuous job\
            \n\t    $ vas_make_cont.py -d {dirname} -j {vjob} -x {partition} -N {nnode} -np {nproc}\
            \n\t    $ vas_make_cont.py -d {dirname} -j {vjob} -x {partition} -N {nnode} -np {nproc} -al : all prepared except POSCAR\
            \n\t::options\
            \n\t    -J for jobname and dirname\
            \n\t    -p for partition: X1-8, X2-12, X3-20 process\
            \n\t    -N number of nodes\
            \n\t    -n number of total processes: {nproc} <= {nnode} * {nXn[partition]}, which proceed -c\
            \n\t    -c (--cpus-per-task: ncpu/2 per node for memory {ncpu}\
            \n\t    --ntasks-per-node {nXn[partition]/2} in case doesnot know ncpu/node\
            "
    slurm.siesta = " === Job submission"
    slurm.siesta += f"\n        sbatch -J {dirname} -p X{partition} -N {nnode} -n {nproc} /home/joonho/sandbox/pypbs/slurm_siesta.sh"
    #slurm.siesta += f"\n        sbatch -J {dirname} -p X{partition} -N {nnode} -n {nproc} --export=exe='gam' /home/joonho/sandbox/pypbs/slurm_sbatch.sh"
    if nproc != nnode * nXn[partition]:
        print("Warning!! Not using all the processes in the node")

    ###### JOB start 
    if work == 'amp':
        ### modify setting here
        settings = {'force':'force',
                    'dtype': 'interval',
                    'dtype-te': 'interval-te',
                    'gs'   : 'log'
                    }
        print(f"{work}:: in {job_submit}: w. {settings['force']} data:{settings['dtype']} gs:{settings['gs']}-scale ")
        if job_submit == 'node':
            print("Train:\namp_run.py -f OUTCAR", amp_node_str[settings['force']], amp_data_type[settings['dtype']], amp_gaussian_function_param[settings['gs']], amp_node_str_suff)
            print("Test:\namp_run.py -f OUTCAR -j te -nc 10", amp_data_type[settings['dtype-te']], amp_node_str_suff)
            print("\n=== New Method")
            print("amp_wrapper.py -js node [-te] [-c| &] ")
            print("\t -c: for check string but not to submit jobs")
            print("GA")
            print("ga_amp_run.py -js node -nch 8 -hl 10 -nn 10 -np 5 [-c] &")
        elif job_submit == 'qsub':
            print("Structure of script to run multiple amp jobs::\
            \n    amp_jobs.sh wrapper_subdir [| sh]; which runs\
            \n    amp_wrapper.py -js qsub -qn Npn${n}hl20 -sf NN${n} -hl 20 20 -dl 1000 1300 3500 3600 &; which runs\
            \n    make mlet_tr(te).csh and run qsub\
            \n    qsub -N HL44tr -pe numa 8 -l mem=12G mlet_tr.csh\
            \n\tmlet_tr.csh should have -f force to make force_derivative directory\
            ")
            print("AMP JOBS")
            print("    amp_jobs.sh job[wrapper wrapper_subdir grep sh etc]")
            print("\nAMP WRAPPER:: train and test with failed-convergence")
            print(f"\tamp_wrapper.py -js qsub [-j [te|tr|trte]] -qn HL44 [-hl 4 4] [-dl int int]  [-s] [-c] &")
            print(f"\tamp_wrapper.py -js qsub &  !!! all options in the script")
            print("\tamp_wrapper.py -js qsub -j te -qn check &")
            print("\tamp_wrapper.py -js qsub -j tr -qn NN9p2 -dl 360 720 &")
            print("\tamp_wrapper.py -js qsub -sf logpn10max200 -qn HL1010 -hl 10 10 -dl 1000 1100 3500 3600 &")
            print("\tqsub -N HL44tr -pe numa 8 -l mem=12G mlet_tr.csh")
            print("    Water molecule:\
                    \n\tamp_wrapper.py -f w1.extxyz -js qsub -qn wat1pn5 -sf NN5 -hl 10 10 -dl 0 300 900 1000 -nt 1000 &\
                    \n\tamp_wrapper.py -f w1.extxyz -qn wat1pn5 -sf NN5 -hl 10 10 -dl 0 700 900 1000 -nt 1000 &\
                    ")

            print("    DB folder")
            print("\tamp_wrapper.py -js qsub -j tr -qn hl4 -dl 0 400 &")
            print("    DB")
            print("\tdiramp.sh fp | sh")
            print("    TRAIN")
            print("\tamp_jobs.sh wrapper")
            print("\tamp_wrapper.py -js qsub -qn HL1010 -hl 10 10 &")
            print("\t    : queue submit job HL44tr-HL44te consecutively")

            print(amp.ga)
        ### this part is for common
        print(amp.clean)
        print(amp.statistics)
        print(amp.analysis)
        #print("\nAnalysis")
        #print("    get_frmse NN 2 | sort -n -t N -k 3")
        if 'kw1' in locals():
            d_pre=kw1
        else:
            d_pre='ND'
        print(plot.first)
        print("    TEST")
        print(f"\tampplot_stat_dir.py -p {d_pre} -le mse -t 'Energy Training' -k hl")
        print("\t    in case dir scan is hl: -k hl")
        print(f"\tampplot_stat_dir.py -p {d_pre} -le mse -t 'Force (Energy-only training)' -yl '(eV/A)^2' Correlation")
        print(f"\tampplot_stat_dir.py -p {d_pre} -le maxres -t 'Force (Energy-only training)' -yl eV/A")
        print(f"\tampplot_stat_dir.py -p {d_pre} -le mse -t 'Energy (Force-Training)'")
        print(f"\tampplot_stat_dir.py -p {d_pre} -le mse -t 'Force (Force-Training)' -yl '(eV/A)^2' Correlation")
        print(f"\tampplot_stat_dir.py -p {d_pre} -le maxres -t 'Force (Force-Training)' -yl eV/A ")
        print("\tampplot_stat_dir.py -p {d_pre} -le maxres -t 'Force (Force-Training)' -yl \"F\$_{rmse}\$ (eV/A)\" \"F\$_{maxres}\$ (eV/A)\"")
        print("\tampplot_stat_dir.py -p {d_pre} -le maxres -xl Nparam -t 'Force (Force-Training)' -yl \"F\$_{rmse}\$ (eV/A)\" \"F\$_{maxres}\$ (eV/A)\"")
        print("\tampplot_stat_dir.py -p NN -le maxres -sd Ndata300 -f f -xl Nparam -t 'Ndata 300 (Gs-pow)' -yl \"F\$_{rmse}\$ (eV/A)\" \"F\$_{maxres}\$ (eV/A)\"")
        print("\tampplot_stat_dir.py -p nd100i -le mse -t 'Energy (ND=100)' -xl \"Data index\"")
        print("\tampplot_stat_dir.py -p nd100i -le maxres -t 'Force (ND=100)' -yl \"F\$_{rmse}\$ (eV/A)\" \"F\$_{maxres}\$ (eV/A)\"  -xl 'Data index'")
    elif work == 'qchem':
        if job_submit == 'qsub':
            print("Run Q-Chem")
            print("    N.B.:")
            print("\tv5.1 parallel runs only at skylake@node!!")
            if not nodename:
                print("    qsub_server.py qchem -qj CC6 -i 6-CC-NiFe.in -n 16 -m 3 [-no skylake@node14]")
            else:
                print(f"    qsub_server.py qchem -qj CC6 -i 6-CC-NiFe.in -n 16 -m 3 -no {nodename}")
            print("\tN.B.: qchem parallel 5.1p doesnot run at Node of (haswell)opt, (sandy)slet -> use skylake@node00")
            print("\t-n for nproc, -m for mem")
            print("NBO analysis")
            print("    qcout_nbo.py nbo -f 5-NiFe.out -a C O O Ni P P Fe P P N -g C O O")

    elif work == 'ga':
        print(ga.amp)
        print(ga.dyn)

    ### run VASP in SLURM
    elif work == 'slurm':
        print(f"Usage for {work}")
        print(f"\tnproc {nproc} = nnode {nnode} * nproc/node_partition {nXn[partition]}")    
        print(f"\t {os.path.basename(__file__)} {work} -j vasp -qn dirname -p {partition} -N {nnode} -n {nproc}")
        print("Run in slurm")
        print("\tsleep:")
        print(f"\t    sbatch -J {qname} -p X{partition} -N {nnode} -n {nproc} --nodelist {nodelist} /home/joonho/sandbox/pypbs/slurm_sbatch_sleep.sh")
        if not subwork:
            print("use -sw subwork [amp|vasp|mldyn|nc|vasnc|vaspnc]")

        elif subwork == 'amp':
            if nhl:
                hlstr1 = list2str(nhl, delimit=" ")
                hl2str = list2str(nhl)
                if not qname:
                    qname = 'hl'+hl2str
            if not qname:
                qname = 'amptest'
            print(f"\tsbatch -J {qname} -p X{partition} -N {nnode} -n {nproc} /home/joonho/sandbox/pypbs/slurm_sbatch_vasp.sh")
            print(f"\tsbatch -J {qname} -p X{partition} -N {nnode} -n {nproc} /home/joonho/sandbox/pypbs/slurm_sbatch_py.sh")
            print(slurm.amp)
            print(f"\t1. amp_jobs.sh fp {inf} {partition} {idata} {ndata} | sh")
            print(f"\t1.9 amp_jobs.sh db {inf} 2 60 | sh")
            print(f"\t2. amp_jobs.sh db {inf} 1 | sh")
            print(f"\t3. amp_jobs.sh wrapper2 {inf} 2 | sh")
        
        elif re.search('vas', subwork) or re.search('nc', subwork):
            if re.search('vas', subwork):
                print(f"=== Usage ===\
                        \n\t{os.path.basename(__file__)} {work} -s {poscar} -p 3 -N 4 -n total_proc\
                        \n\tN.B.:: (POSCAR.)name is dirname and qname in qsub") 
                print("=== Prepare VASP directory in platinum")
                print(f"    (1) vas_make_ini.py -s {poscar}")
                print("\tmake POSCAR POTCAR KPOINTS INCAR & directory")
                print(f"    python -m myvasp -j getmag -p {poscar}")
                print(f"    sed -i 's/.*MAGMOM.*/ mag_moment/' INCAR")
                print("\tto modify MAGMON in INCAR from POSCAR in module myvasp.py")
                print(slurm.vas)
            if re.search('nc', subwork):
                print("OORINano/catalysis with SBATCH:")
                if keyvalues:
                    kv  = keyvalues[0]
                else:
                    kv  =   'orr'
                
                #if not qname:
                #    qname=f'{keyvalues}_test'
                print(f"\t:: cat=['orr' (default),'her'], pos=['cp' (default: copy), 'any_other_word']")
                print(f"\tsbatch -J {qname} -p X{partition} -N {nnode} -n {nproc}  slurm_sbatch_nc.sh")
                print(f"\tsbatch -J {qname} -p X{partition} -N {nnode} -n {nproc} --export=cat='{kv}' --export=pos='cp' slurm_sbatch_nc.sh")
                print(f"\t:: if not 'cp', generate")
                print(f"\tsbatch -J {qname} -p X{partition} -N {nnode} -n {nproc} --export=pos='gen' slurm_sbatch_nc.sh")
                print("OORINano/qtnegf with SBATCH:")
                print(f"\t:: Make subdirectory with {qname} - default")
                print(f"\tsbatch -J {qname} -p X{partition} -N {nnode} -n {nproc}  slm_qtnegf.sh")
                print(f"\t:: Run in main directory")
                print(f"\tsbatch -J {qname} -p X{partition} -N {nnode} -n {nproc}  --export=sub=0 slm_qtnegf.sh")
                print("\nOORINano Package Development:")
                print(nc.build)
        elif subwork == 'mldyn':
            print("SBATCH:")
            print(f"\tsbatch -J {qname} -p X2 -N 1 -n 1 --export=hl='{hlstr1}',sf='hl{hl2str}.pt' /home/joonho/sandbox/pypbs/slurm_sbatch_py.sh")
            print("\t    slurm_sbatch_py.sh")
            print("\t\tml_lorenz.py tr -hl $hl -ms $sf")
            print("Direct run:")
            print("    GPU is coded in PyTorch")
            print("\tssh pt|fe; ssh n076 (GPU). fe also has GPU")
            print("    mllorenz.py --> ml_lorenz_gpu.py")
            print("    Loop for nhls")
            print("\tgpu_ml_batch.sh")
            print(f"    mllorenz.py tr -hl {hlstr1} -ms {qname}.pt")
            print(f"    mllorenz.py te -m {qname}.pt -dbp 2")
            print(f"    mllorenz.py -> ml_lorenz_gpu.py|ml_lorenz_cpu.py")
            print("\t-hl hidden layers")
            print("\t-m load saved model of a.pt")
            print("\t-ms save pytorch model such as hl101010")
            print("\t-mi max iteration, default=10^5")
            print("\t-dbp partition [2|3]: data to tr, [val, and] te")
            print(ga.dyn)
        elif subwork == 'crr':
            print("    CRR:DAC =========================")
            print(nc.build)
            print(dac.build)
        print(slurm.pbs)
    elif work == 'kisti' or work == 'pbs':
        print("=== PBS in KISTI ===")
        print("    Bash ")
        print(kisti.shell)
        print(kisti.pbs)
        print("=== VASP in KISTI ===")
        if subwork == 'vasp':
            print(comment_subj.vasp.run)
            print(comment_sys.server.kisti.pbs)
        print(kisti.vas)
        print("=== NanoCore in KISTI ===")
        print("qsub")
        if keyvalues:
            kv  = keyvalues[0]
        else:
            kv  =   'orr'
        
        if not qname:
            qname=f'{keyvalues}_test'
        print(f"\tqsub -N {qname} pbs_vasp_kisti_skl.sh")
        print(f"\tqsub -N {qname} -v cat='{kv}' pbs_vasp_kisti_skl.sh")
        print("\nNonoCore Package Development:")
        print(nc.build)

    elif work == 'vasp':
        print("KISTI: vasp")
        print(kisti.vas)
        print(kisti.shell)
        print(kisti.pbs)
        print("Pt(platinum:slurm): vasp")
        print(slurm.vas)
    elif work == 'siesta':
        print("Pt(platinum:slurm): siesta")
        print(slurm.siesta)
    else:
        print("build more works")
    return 0 
    

def main():

    parser = argparse.ArgumentParser(description="show command Amp/Qchem/ etc ")
    parser.add_argument('work', choices=['slurm','kisti','amp','qchem','ga','gpu','vasp','pbs','siesta'],  help="one of amp, qchem, mldyn for ML dyn")
    parser.add_argument('-sw', '--subwork', choices=['vasp', 'mldyn', 'nc', 'crr', 'amp', 'vaspnc', 'vasnc'], help="one of amp, qchem, mldyn for ML dyn")
    parser.add_argument('-js','--job_submit', default='qsub', choices=['chi','qsub','getqsub', 'node'],  help="where the job running ")
    parser.add_argument('-qn', '-q', '--qname', default='test', help="queue name for qsub shown by qstat")
    parser.add_argument('-vj','--vas_job', default='sp', help="vasp job input for vas_make_ini.py")
    parser.add_argument('-sj','--vas_subjob', default='mag', help="vasp subjob input for vas_make_ini.py to make INCAR.spmag")
    parser.add_argument('-kv', '--keyvalues', nargs='*', help='change a keyword in print')
    parser.add_argument('-no', '--nodename', help='if needed, specify nodename')
    ### flowing slurm option
    parser.add_argument('-inf', '--infile', help='input file in case')
    parser.add_argument('-N', '--nnode', default=1, type=int, help='number of nodes: if needed')
    parser.add_argument('-np', '--nproc', type=int, help='number of process: if needed')
    parser.add_argument('-nl', '--nodelist',  help='node list to assign nodes')
    parser.add_argument('-x', '--xpartition', default=3, type=int, choices=[1,2,3,4,5,6], help='if needed, specify nodename')
    parser.add_argument('-s', '--poscar', default='POSCAR.name', help='if needed, specify nodename')
    parser.add_argument('-id', '--idata', default=0, type=int, help='start index of data')
    parser.add_argument('-nd', '--ndata', default=100, type=int, help='amount of data')
    mlg = parser.add_argument_group(title = 'machine learning args')
    mlg.add_argument('-dt', '--data_type', default='int', choices=['int','div'], help="data selection type")
    mlg.add_argument('-ft', '--func_type', default='log', choices=['log','pow'], help="gaussian symmetry function parameter type")
    mlg.add_argument('-hl', '--hidden_layers', nargs='*', default=['4','4','4'], help="hidden layers in integer")

    args = parser.parse_args()

    if not args.infile:
        if args.subwork == 'amp':
            infile = 'OUTCAR'
        elif args.subwork == 'nc':
            infile = 'test_HER.py'
        else:
            infile = 'OUTCAR'
    else:
        infile = args.infile

    if 'POSCAR' in args.qname:
        qname=args.qname[7:]
    else:
        qname=args.qname

    show_command(args.work,args.subwork,args.job_submit,qname,args.vas_job,args.vas_subjob,infile,args.keyvalues,args.nodename,args.nnode,args.nproc,args.nodelist, args.func_type,args.data_type,args.xpartition,args.poscar, args.hidden_layers, args.idata, args.ndata)

if __name__ == "__main__":
    main()
