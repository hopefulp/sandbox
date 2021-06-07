#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_files, MyClass, list2str

amp     =   MyClass('amp')
plot    =   MyClass('plot')
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
amp_qsub_str_suff=" /home/joonho/sandboxg/pypbs/sge_amp.csh"                        

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
plot.first      =   "\nPLOT\
                    \n    TRAIN & TEST\
                    \n\tDraw file\
                    \n\t    f2d_plot.py train.dat -pd -xs 2\
                    \n\t    f2d_plot.py train.dat test.dat -pd -t 'TR-100:TE-100'\
                    \n\t\t-pd pandas, -xs xspacing, -xst xspacing type (default:numerically)\
                    \n\tampdir.sh grep\
                    \n\tampplot_dir.py -p NN -t 'Ndata 100 (Gs-pow)'\
                    \n\tampplot_dir.py -p NN -t 'Training Set: Gs-pow' -y tr -yd . Nd300hl2020 Nd500hl2020\
                    \n\tampplot_dir.py -p NN -t 'Test Set    : Gs-pow' -y te -yd . Nd300hl2020 Nd500hl2020\
                    "

def show_command(job, subjob, job_submit, qname, keyvalues, nodename, nnode, nproc, sftype, dtype, partition,poscar, nhl):
    
    if keyvalues:
        kw1 = keyvalues[0]

    if job == 'amp':
        ### modify setting here
        settings = {'force':'force',
                    'dtype': 'interval',
                    'dtype-te': 'interval-te',
                    'gs'   : 'log'
                    }
        print(f"{job}:: in {job_submit}: w. {settings['force']} data:{settings['dtype']} gs:{settings['gs']}-scale ")
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
            \n    ampdir.sh wrapper_subdir [| sh]; which runs\
            \n    amp_wrapper.py -js qsub -qn Npn${n}hl20 -sf NN${n} -hl 20 20 -dl 1000 1300 3500 3600 &; which runs\
            \n    make mlet_tr(te).csh and run qsub\
            \n    qsub -N HL44tr -pe numa 8 -l mem=12G mlet_tr.csh\
            \n\tmlet_tr.csh should have -f force to make force_derivative directory\
            ")
            print("AMP JOBS")
            print("    ampdir.sh job[wrapper wrapper_subdir grep sh etc]")
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
            print("\tampdir.sh wrapper")
            print("\tamp_wrapper.py -js qsub -qn HL1010 -hl 10 10 &")
            print("\t    : queue submit job HL44tr-HL44te consecutively")

            print("GA")
            print("ga_amp_run.py -js qsub -nch 12 -hl 7 -nn 15 -np 4")
        elif job_submit == 'getqsub':
            print(f"sge_amp.py -qj {qname} -m 3G", amp_getqsub_str[settings['force']], amp_getqsub_datatype[settings['dtype']], amp_getqsub_symmfuncttype[settings['gs']]) 
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
    elif job == 'qchem':
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

    ### run VASP in SLURM
    elif job == 'slurm':
        print("Run in slurm")
        if nhl:
            hlstr1 = list2str(nhl, delimit=" ")
            hl2str = list2str(nhl)
            if not qname:
                qname = 'hl'+hl2str
        if not qname:
            qname = 'amptest'
        print("\tsbatch -J {qname} -p X1 -N 4 -n 96 /home/joonho/sandbox_gl/pypbs/slurm_sbatch_vasp.sh")
        print(f"\tsbatch -J {qname} -p X2 -N 1 -n 1 /home/joonho/sandbox_gl/pypbs/slurm_sbatch_py.sh")
        if not subjob:
            print("use -sj subjob [vasp|mldyn]")
        elif subjob == 'vasp':
            print(f"=== Usage ===\
                    \n\t{os.path.basename(__file__)} {job} -s {poscar} -p 3 -N 4 -n total_proc\
                    \n\tN.B.:: (POSCAR.)name is dirname and qname in qsub") 
            print("=== Prepare VASP directory in platinum")
            print(f"    (1) vas_make_ini.py -s {poscar}")
            print("\tmake POSCAR POTCAR KPOINTS INCAR & directory")
            #if poscar:
            #    com = f"python -m myvasp -j getmag -p {poscar}"
            print(f"    python -m myvasp -j getmag -p {poscar}")
            #mag_moment = os.system(com)
            print(f"    sed -i 's/.*MAGMOM.*/ mag_moment/' INCAR")
            print("\tto modify MAGMON in INCAR from POSCAR in module myvasp.py")
            if poscar and re.match('POSCAR', poscar) :
                if poscar[7:] != 'name':
                    dirname = poscar[7:]
            if 'dirname' not in locals():
                dirname = qname
            nXn = { '1': 8, '2': 12, '3': 20, '4':24, '5':32 }
            if nnode == '1':
                nproc = nXn[partition]
                print(f"nproc = {nproc}")
            print("=== Job submission")
            print(f"    (2) sbatch -J {dirname} -p X{partition} -N {nnode} -n {nproc} /home/joonho/sandbox_gl/pypbs/slurm_sbatch.sh")
            print("\toptions::")
            print("\t    -J for jobname and dirname")
            print("\t    -p for partition: X1-8, X2-12, X3-20 process")
            print("\t    -N number of nodes ")
            print(f"\t    -n number of total processes: {nproc} <= {nnode} * {nXn[partition]}")
            if nproc != nnode * nXn[partition]:
                print("Warning!! Not using all the processes in the node")
            print("=== Node check ===")
            print("\tpe     to get free processes")
            print("\tpef    to get free nodes")
            print("\tpestat to see all nodes")
            print("\tqstat -u joonho     to check my job")
        elif subjob == 'mldyn':
            print("SBATCH:")
            print(f"\tsbatch -J {qname} -p X2 -N 1 -n 1 --export=hl='{hlstr1}',sf='hl{hl2str}.pt' /home/joonho/sandbox_gl/pypbs/slurm_sbatch_py.sh")
            print("\t    slurm_sbatch_py.sh")
            print("\t\tml_lorenz.py tr -hl $hl -ms $sf")
            print("Direct run:")
            print(f"    ml_lorenz.py tr -hl {hlstr1} -ms {qname}.pt")
            print(f"    ml_lorenz.py te -m {qname}.pt -dbp 2")
            print("\t-hl hidden layers")
            print("\t-m load saved model of a.pt")
            print("\t-ms save pytorch model such as hl101010")
            print("\t-mi max iteration, default=10^5")
            print("\t-dbp partition [2|3]: data to tr, [val, and] te")

    else:
        print("build more jobs")
    return 0 
    

def main():

    parser = argparse.ArgumentParser(description="show command Amp/Qchem/ etc ")
    parser.add_argument('job', choices=['amp','qchem','slurm'],  help="one of amp, qchem, mldyn for ML dyn")
    parser.add_argument('-sj', '--sub_job', choices=['vasp', 'mldyn'],  help="one of amp, qchem, mldyn for ML dyn")
    parser.add_argument('-js','--job_submit', default='qsub', choices=['chi','qsub','getqsub', 'node'],  help="where the job running ")
    parser.add_argument('-qn', '--qname', help="queue name for qsub shown by qstat")
    parser.add_argument('-k', '--keyvalues', nargs='*', help='change a keyword in print')
    parser.add_argument('-no', '--nodename', help='if needed, specify nodename')
    ### flowing slurm option
    parser.add_argument('-N', '--nnode', default=1, type=int, help='number of nodes: if needed')
    parser.add_argument('-n', '--nproc', default=1, type=int, help='number of process: if needed')
    parser.add_argument('-p', '--partition', default='1', choices=['1','2','3','4','5'], help='if needed, specify nodename')
    parser.add_argument('-s', '--poscar', default='POSCAR.name', help='if needed, specify nodename')
    mlg = parser.add_argument_group(title = 'machine learning args')
    mlg.add_argument('-dt', '--data_type', default='int', choices=['int','div'], help="data selection type")
    mlg.add_argument('-ft', '--func_type', default='log', choices=['log','pow'], help="gaussian symmetry function parameter type")
    mlg.add_argument('-hl', '--hidden_layers', nargs='*', help="hidden layers in integer")

    args = parser.parse_args()

    show_command(args.job,args.sub_job,args.job_submit,args.qname,args.keyvalues, args.nodename,args.nnode,args.nproc, args.func_type,args.data_type,args.partition,args.poscar, args.hidden_layers)

if __name__ == "__main__":
    main()
